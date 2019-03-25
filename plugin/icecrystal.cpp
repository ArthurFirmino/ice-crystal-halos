/*
	This file is part of ice-crystal-halos, a rendering plugin for Mitsuba.

	Copyright (c) 2019 by Arthur Firmino.
*/

#include <mitsuba/core/frame.h>
#include <mitsuba/render/phase.h>
#include <mitsuba/render/medium.h>
#include <mitsuba/render/sampler.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/warp.h>

#include <array>

MTS_NAMESPACE_BEGIN

#define ICECRYSTAL_DEBUG 0

class IcecrystalPhaseFunction : public PhaseFunction {
public:

	IcecrystalPhaseFunction(const Properties &props) : PhaseFunction(props) {

		// Open phase function file
		FileResolver *fResolver = Thread::getThread()->getFileResolver();
		fs::path path = fResolver->resolve(props.getString("filename"));
		if (!fs::exists(path))
			Log(EError, "Icecrystal phase function file \"%s\" could not be found!",
				path.string().c_str());
		fs::ifstream in(path, std::ios::binary);

		// Verify number of spectral samples matches
		uint32_t spectralSamples;
		in.read(reinterpret_cast<char*>(&spectralSamples), sizeof(uint32_t));
		if (spectralSamples == 1) {
			m_monochromatic = true;
		} else if (spectralSamples == SPECTRUM_SAMPLES) {
			m_monochromatic = false;
		} else {
			Log(EError, "Mismatched number of phase function samples!");
		}

		// Read in table dimensions
		in.read(reinterpret_cast<char*>(&m_nThetaI), sizeof(uint32_t));
		in.read(reinterpret_cast<char*>(&m_nThetaO), sizeof(uint32_t));
		in.read(reinterpret_cast<char*>(&m_nDeltaPhi), sizeof(uint32_t));

		// Read in table
		m_phaseFunctionArray.resize(m_nThetaI*m_nThetaO*m_nDeltaPhi);
		in.read(reinterpret_cast<char*>(&m_phaseFunctionArray[0]), m_phaseFunctionArray.size() * sizeof(Float));
		m_sigmaScaleArray.resize(m_nThetaI);
		in.read(reinterpret_cast<char*>(&m_sigmaScaleArray[0]), m_sigmaScaleArray.size() * sizeof(Float));

		// Read in spectral phase function if applicable
		if (!m_monochromatic) {
			m_spectralPhaseFunctionArray.resize(m_phaseFunctionArray.size());
			in.read(reinterpret_cast<char*>(&m_spectralPhaseFunctionArray[0]), 
				m_spectralPhaseFunctionArray.size() * SPECTRUM_SAMPLES * sizeof(Float));
		}

		configure();
	}

	IcecrystalPhaseFunction(Stream *stream, InstanceManager *manager)
		: PhaseFunction(stream, manager) {

		m_monochromatic = stream->readBool();
		m_nThetaI = stream->readUInt();
		m_nThetaO = stream->readUInt();
		m_nDeltaPhi = stream->readUInt();

		m_phaseFunctionArray.resize(m_nThetaI*m_nThetaO*m_nDeltaPhi);
		stream->readFloatArray(&m_phaseFunctionArray[0], m_phaseFunctionArray.size());
		m_sigmaScaleArray.resize(m_nThetaI);
		stream->readFloatArray(&m_sigmaScaleArray[0], m_sigmaScaleArray.size());

		if (!m_monochromatic) {
			m_spectralPhaseFunctionArray.resize(m_phaseFunctionArray.size());
			stream->readFloatArray(reinterpret_cast<Float*>(&m_spectralPhaseFunctionArray[0]),
				m_spectralPhaseFunctionArray.size() * SPECTRUM_SAMPLES);
		}

		configure();
	}

	virtual ~IcecrystalPhaseFunction() { }

	void configure() {

		PhaseFunction::configure();

		m_type = ENonSymmetric | EAnisotropic;
		
		// Compute Marginal and Cumulative Distribution Functionns for sampling
		m_lowNThetaI = 32;
		m_lowNThetaO = 32;
		m_lowNDeltaPhi = 32;
		SAssert(!(m_nThetaI % m_lowNThetaI) && !(m_nThetaO % m_lowNThetaO) && !(m_nDeltaPhi % m_lowNDeltaPhi));

		m_marginalCDF.resize(m_lowNThetaI * m_lowNThetaO, 0);
		m_cumulativeDF.resize(m_lowNThetaI * m_lowNThetaO * m_lowNDeltaPhi, 0);

		auto mcdf = [&](int i, int j) { return &m_marginalCDF[i*m_lowNThetaO + j ]; };
		auto cdf = [&](int i, int j, int k) { return &m_cumulativeDF[i*(m_lowNThetaO * m_lowNDeltaPhi) + j*m_lowNDeltaPhi + k]; };

		int nI = m_nThetaI / m_lowNThetaI;
		int nJ = m_nThetaO / m_lowNThetaO;
		int nK = m_nDeltaPhi / m_lowNDeltaPhi;

		auto sumblock = [&](int i, int j, int k) {
			Float I = 0;
			Float W = 0;
			for (int ii = 0; ii < nI; ++ii) {
				for (int jj = 0; jj < nJ; ++jj) {
					for (int kk = 0; kk < nK; ++kk) {
						Float w;
						if (i + ii == 0 || j + jj == 0 || k + kk == 0 || i + ii == m_nThetaI - 1 || j + jj == m_nThetaO - 1 || k + kk == m_nDeltaPhi) w += 0.5f;
						else w += 1.0f;
						I += w*pf(i + ii, j + jj, k + kk);
						W += w;
					}
				}
			}
			return I / W;
		};

		for (unsigned i = 0; i < m_lowNThetaI; ++i) {
			for (unsigned j = 0; j < m_lowNThetaO; ++j) {
				for (unsigned k = 0; k < m_lowNDeltaPhi; ++k) {
					*cdf(i, j, k) += sumblock(i*nI,j*nJ,k*nK);
					if (k > 0) *cdf(i,j,k) += *cdf(i, j, k - 1);
				}
				*mcdf(i, j) = *cdf(i, j, m_lowNDeltaPhi - 1);
				for (unsigned k = 0; k < m_lowNDeltaPhi; ++k) *cdf(i, j, k) /= *cdf(i, j, m_lowNDeltaPhi-1);
				if(j > 0) *mcdf(i, j) += *mcdf(i, j - 1);
			}
			for (unsigned j = 0; j < m_lowNThetaO; ++j) *mcdf(i, j) /= *mcdf(i, m_lowNThetaO - 1);
		}
	}

	void serialize(Stream *stream, InstanceManager *manager) const {

		PhaseFunction::serialize(stream, manager);
		stream->writeBool(m_monochromatic);
		stream->writeUInt(m_nThetaI);
		stream->writeUInt(m_nThetaO);
		stream->writeUInt(m_nDeltaPhi);

		stream->writeFloatArray(&m_phaseFunctionArray[0], m_phaseFunctionArray.size());
		stream->writeFloatArray(&m_sigmaScaleArray[0], m_sigmaScaleArray.size());

		if (!m_monochromatic) {
			stream->writeFloatArray(reinterpret_cast<const Float*>(&m_spectralPhaseFunctionArray[0]),
				m_spectralPhaseFunctionArray.size() * SPECTRUM_SAMPLES);
		}
	}

	Float eval(const PhaseFunctionSamplingRecord &pRec) const {
		
#if ICECRYSTAL_DEBUG
		if (pRec.mRec.orientation.isZero()) std::abort();
		SAssert(!pRec.mRec.orientation.isZero());
#endif
		Frame frame(pRec.mRec.orientation);
		Vector wi = frame.toLocal(pRec.wi);
		Vector wo = frame.toLocal(pRec.wo);
		wi = normalize(wi);
		wo = normalize(wo);

		int idx[3];
		Float T[3];
		getGridCoordinates(wi, wo, idx, T);
		Float ret = 0;
		
		Float b[3][4];
		for (int i = 0; i < 3; ++i) {
			Float t = T[i], t2 = t * t, t3 = t2 * t;
			b[i][0] = -0.5f*t + t2 - 0.5f*t3;
			b[i][1] = 1.0f - 2.5f*t2 + 1.5f*t3;
			b[i][2] = 0.5f*t + 2.0f*t2 - 1.5f*t3;
			b[i][3] = -0.5f*t2 + 0.5f*t3;
		}

		for (int i = -1; i < 3; ++i) {
			Float B = b[0][i + 1];
			for (int j = -1; j < 3; ++j) {
				Float temp = B;
				B *= b[1][j + 1];
				for (int k = -1; k < 3; ++k) {
					ret += (B * b[2][k + 1]) * pf(idx[0] + i, idx[1] + j, idx[2] + k);
				}
				B = temp;
			}
		}

		return ret > 0 ? ret : 0;
	}

	Float sample(PhaseFunctionSamplingRecord &pRec, Sampler *sampler) const {
#if ICECRYSTAL_DEBUG
		SAssert(!pRec.mRec.orientation.isZero());
#endif
		Frame frame(pRec.mRec.orientation);
		Vector wi = frame.toLocal(pRec.wi);

		auto bisearch = [&](const std::function<Float(int)>& f, Float x, int a, int b, Float& negoff) {
			while (a < b) {
				int m = (a + b) / 2;
				if (x <= f(m)) {
					b = m;
				} else {
					if (a == m) break;
					a = m;
				}
			}
			Float f0 = b > 0 ? f(b - 1) : 0;
			Float f1 = f(b);
			negoff = (x - f1) / (f1 - f0);
			return b;
		};

		auto s = sampler->next2D();
		int i = static_cast<int>(m_lowNThetaI * (wi.z + 1.0f) * 0.5f);
		i = i < 0 ? 0 : i >= static_cast<int>(m_lowNThetaI) ? m_lowNThetaI - 1 : i;

		auto mcdf = [&](int j) { return m_marginalCDF[i*m_lowNThetaO + j]; };
		Float off;
		int j = bisearch(mcdf, s.x, 0, m_lowNThetaO - 1, off);
		Float r = 2 * ((1 + j + off) / m_lowNThetaO) - 1.0f;
		Float theta = std::acos( r );

		auto cdf = [&](int k) { return m_cumulativeDF[i*(m_lowNThetaO * m_lowNDeltaPhi) + j * m_lowNDeltaPhi + k]; };
		int k = bisearch(cdf, s.y, 0, m_lowNDeltaPhi - 1, off);
		Float dPhi = M_PI * (k + off) / (m_lowNDeltaPhi - 1.0f);
		if (sampler->next1D() < 0.5) dPhi *= -1;

		auto SphericalPhi = [](const Vector &v) {
			Float p = std::atan2(v.y, v.x);
			return (p < 0) ? (p + 2 * M_PI) : p;
		};
		Float phi = SphericalPhi(wi) + dPhi;

		pRec.wo = frame.toWorld(Vector(cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta)));
		return 1.f;
	}

	Float sample(PhaseFunctionSamplingRecord &pRec, Float &pdf, Sampler *sampler) const {
		if (sample(pRec, sampler) == 0) {
			pdf = 0.f; return 0.f;
		}
		pdf = eval(pRec);
		return 1.f;
	}

	Float pdf(const PhaseFunctionSamplingRecord &pRec) const {
		return eval(pRec);
	}

	Spectrum spectralEval(const PhaseFunctionSamplingRecord &pRec) const {
#if ICECRYSTAL_DEBUG
		SAssert(!pRec.mRec.orientation.isZero());
		SAssert(!m_monochromatic);
#endif

		Frame frame(pRec.mRec.orientation);
		Vector wi = frame.toLocal(pRec.wi);
		Vector wo = frame.toLocal(pRec.wo);

		int idx[3];
		Float T[3];
		getGridCoordinates(wi, wo, idx, T);
		Spectrum ret(0.0f);

		Float b[3][4];
		for (int i = 0; i < 3; ++i) {
			Float t = T[i], t2 = t * t, t3 = t2 * t;
			b[i][0] = -0.5f*t + t2 - 0.5f*t3;
			b[i][1] = 1.0f - 2.5f*t2 + 1.5f*t3;
			b[i][2] = 0.5f*t + 2.0f*t2 - 1.5f*t3;
			b[i][3] = -0.5f*t2 + 0.5f*t3;
		}

		for (int i = -1; i < 3; ++i) {
			Float B = b[0][i+1];
			for (int j = -1; j < 3; ++j) {
				Float temp = B;
				B *= b[1][j+1];
				for (int k = -1; k < 3; ++k) {
					ret += (B * b[2][k+1]) * spectralPf(idx[0] + i, idx[1] + j, idx[2] + k);
				}
				B = temp;
			}
		}
		ret.clampNegative();
#if ICECRYSTAL_DEBUG
		if (!ret.isValid()) {
			std::string s = ret.toString();
			Log(EError, "Invalid interpolated spectrum: %s", ret.toString().c_str());
		}
#endif
		return ret;
	}

	Spectrum spectralSample(PhaseFunctionSamplingRecord &pRec,
		Sampler *sampler) const {
		Float pdf;
		return spectralSample(pRec, pdf, sampler);
	}

	Spectrum spectralSample(PhaseFunctionSamplingRecord &pRec,
		Float &pdf, Sampler *sampler) const {
#if ICECRYSTAL_DEBUG
		SAssert(!m_monochromatic);
#endif
		sample(pRec, pdf, sampler);
		if (pdf <= 0) {
			pdf = 0;
			return Spectrum(0.f);
		}
		return spectralEval(pRec) / pdf;
	}

	bool needsDirectionallyVaryingCoefficients() const { return true; }

	Float sigmaDir(Float cosTheta) const {

		cosTheta = std::fabsf(cosTheta);
		int i = std::min(int((m_nThetaI - 1) * cosTheta), static_cast<int>(m_nThetaI) - 2);
		float t = (m_nThetaI - 1) * cosTheta - i, t2 = t * t, t3 = t2 * t;

		Float b[4];
		b[0] = -0.5f*t + t2 - 0.5f*t3;
		b[1] = 1.0f - 2.5f*t2 + 1.5f*t3;
		b[2] = 0.5f*t + 2.0f*t2 - 1.5f*t3;
		b[3] = -0.5f*t2 + 0.5f*t3;

		int i0 = i == 0 ? i + 1 : i - 1;
		int i3 = i == (static_cast<int>(m_nThetaI) - 2) ? i : i + 2;

		return b[0] * m_sigmaScaleArray[i0] + b[1] * m_sigmaScaleArray[i] + b[2] * m_sigmaScaleArray[i + 1] + b[3] * m_sigmaScaleArray[i3];
	}

	Float sigmaDirMax() const {
		return 1.f;
	}

	bool isSpectral() const {
		return !m_monochromatic;
	}

	std::string toString() const {
		return "IcecrystalPhaseFunction";
	}

	MTS_DECLARE_CLASS()
private:

	bool m_monochromatic;
	uint32_t m_nThetaI, m_nThetaO, m_nDeltaPhi;
	
	std::vector<Float> m_phaseFunctionArray;
	std::vector<Float> m_sigmaScaleArray;
	std::vector<std::array<Float, SPECTRUM_SAMPLES>> m_spectralPhaseFunctionArray;

	uint32_t m_lowNThetaI, m_lowNThetaO, m_lowNDeltaPhi;
	std::vector<Float> m_marginalCDF;
	std::vector<Float> m_cumulativeDF;

	void getGridCoordinates(const Vector &localWi, const Vector &localWo, int *idx, float *T) const {
		
		auto SphericalPhi = [](const Vector &v) {
			Float p = std::atan2(v.y, v.x);
			return (p < 0) ? (p + 2 * M_PI) : p;
		};

		auto AngleDifference = [](Float a, Float b) {
			Float r = fmod(b - a, 2.f*M_PI);
			if (r < -M_PI)
				r += 2.f*M_PI;
			if (r >= M_PI)
				r -= 2.f*M_PI;
			return r;
		};

		Float X[3];
		X[0] = localWi.z;
		X[1] = localWo.z;
		X[2] = abs(AngleDifference(SphericalPhi(localWi), SphericalPhi(localWo)));
		
		Float bounds[3][2] = { {-1.f,1.f}, {-1.f,1.f}, {0,M_PI} };
		int N[3]; N[0] = m_nThetaI; N[1] = m_nThetaO; N[2] = m_nDeltaPhi;

		for (int i = 0; i < 3; ++i) {
#if ICECRYSTAL_DEBUG
			if (!(X[i] >= bounds[i][0] && X[i] <= bounds[i][1])) {
				Log(EWarn, "x[%i] = %f out of bounds [ %f , %f ]", i, X[i], bounds[i][0], bounds[i][1]);
				std::abort();
			}
#endif
			Float x = (X[i] - bounds[i][0]) / (bounds[i][1] - bounds[i][0]);
			idx[i] = std::min(int(x * (N[i] - 1)), N[i] - 2);
			T[i] = (N[i] - 1)*x - idx[i];
#if ICECRYSTAL_DEBUG
			if (!(T[i] >= 0 && T[i] <= 1)) {
				Log(EError, "T = %f outside of [0,1], i = %i, X = %f, idx = %i, x = %f",
					T[i], i, X[i], idx[i], x);
			}
#endif
		}
	}

	Float pf(int i, int j, int k) const {
		if (i < 0) i = -i;
		if (j < 0) j = -j;
		if (k < 0) k = -k;
		if (i > static_cast<int>(m_nThetaI) - 1) i = 2 * (m_nThetaI - 1) - i;
		if (j > static_cast<int>(m_nThetaO) - 1) j = 2 * (m_nThetaO - 1) - j;
		if (k > static_cast<int>(m_nDeltaPhi) - 1) k = 2 * (m_nDeltaPhi - 1) - k;
		return m_phaseFunctionArray[ i*(m_nThetaO * m_nDeltaPhi) + j*m_nDeltaPhi + k ];
	}

	Spectrum spectralPf(int i, int j, int k) const {
		if (i < 0) i = -i;
		if (j < 0) j = -j;
		if (k < 0) k = -k;
		if (i > static_cast<int>(m_nThetaI) - 1) i = 2 * (m_nThetaI - 1) - i;
		if (j > static_cast<int>(m_nThetaO) - 1) j = 2 * (m_nThetaO - 1) - j;
		if (k > static_cast<int>(m_nDeltaPhi) - 1) k = 2 * (m_nDeltaPhi - 1) - k;
		return Spectrum(const_cast<Float*>(reinterpret_cast<const Float*>(&m_spectralPhaseFunctionArray[i*(m_nThetaO * m_nDeltaPhi) + j * m_nDeltaPhi + k])));
	}
};

MTS_IMPLEMENT_CLASS_S(IcecrystalPhaseFunction, false, PhaseFunction)
MTS_EXPORT_PLUGIN(IcecrystalPhaseFunction, "Icecrystal phase function");
MTS_NAMESPACE_END