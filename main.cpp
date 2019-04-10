/*
 *  Copyright (c) 2019 Dmitry Kazakov <dimula73@gmail.com>
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <QCoreApplication>

#ifdef _WIN64
#define MEMALIGN_ALLOC(p, a, s) ((*(p)) = _aligned_malloc((s), (a)), *(p) ? 0 : errno)
#define MEMALIGN_FREE(p) _aligned_free((p))
#else
#define MEMALIGN_ALLOC(p, a, s) ((*(p)) = ::aligned_alloc((a), (s)), *(p) ? 0 : errno)
#define MEMALIGN_FREE(p) std::free((p))
#endif

#include "KoStreamedMath.h"
#include <QElapsedTimer>
#include <QDebug>

struct OptiDiv {
    static ALWAYS_INLINE float divScalar(const float& divident, const float& divisor) {
#ifdef __SSE__
        float result;

        __m128 x = _mm_set_ss(divisor);
        __m128 y = _mm_set_ss(divident);
        x = _mm_rcp_ss(x);
        x = _mm_mul_ss(x, y);


        _mm_store_ss(&result, x);
        return result;
#else
        return divident / divisor;
#endif

    }

    static ALWAYS_INLINE float_v divVector(float_v::AsArg divident, float_v::AsArg  divisor) {
#ifdef __SSE__
        return divident * Vc::reciprocal(divisor);
#else
        return divident / divisor;
#endif

    }
};


// \see docs in AlphaDarkenCompositor32
template <bool haveMask, bool src_aligned>
ALWAYS_INLINE void compositeVectorOver(const quint8 *src, quint8 *dst, const quint8 *mask, float opacity)
{
    float_v src_alpha;
    float_v dst_alpha;

    src_alpha = KoStreamedMath::template fetch_alpha_32<src_aligned>(src);

    bool haveOpacity = opacity != 1.0;
    float_v opacity_norm_vec(opacity);

    float_v uint8Max((float)255.0);
    float_v uint8MaxRec1((float)1.0 / 255.0);
    float_v zeroValue(Vc::Zero);
    float_v oneValue(Vc::One);

    src_alpha *= opacity_norm_vec;

    if (haveMask) {
        float_v mask_vec = KoStreamedMath::fetch_mask_8(mask);
        src_alpha *= mask_vec * uint8MaxRec1;
    }

    // The source cannot change the colors in the destination,
    // since its fully transparent
    if ((src_alpha == zeroValue).isFull()) {
        return;
    }

    dst_alpha = KoStreamedMath::template fetch_alpha_32<true>(dst);

    float_v src_c1;
    float_v src_c2;
    float_v src_c3;

    float_v dst_c1;
    float_v dst_c2;
    float_v dst_c3;


    KoStreamedMath::template fetch_colors_32<src_aligned>(src, src_c1, src_c2, src_c3);
    float_v src_blend;
    float_v new_alpha;

    if ((dst_alpha == uint8Max).isFull()) {
        new_alpha = dst_alpha;
        src_blend = src_alpha * uint8MaxRec1;
    } else if ((dst_alpha == zeroValue).isFull()) {
        new_alpha = src_alpha;
        src_blend = oneValue;
    } else {
        /**
         * The value of new_alpha can have *some* zero values,
         * which will result in NaN values while division. But
         * when converted to integers these NaN values will
         * be converted to zeroes, which is exactly what we need
         */
        new_alpha = dst_alpha + (uint8Max - dst_alpha) * src_alpha * uint8MaxRec1;

        // Optimized version of:
        //     src_blend = src_alpha / new_alpha;
        src_blend = OptiDiv::divVector(src_alpha, new_alpha);

    }

    if (!(src_blend == oneValue).isFull()) {
        KoStreamedMath::template fetch_colors_32<true>(dst, dst_c1, dst_c2, dst_c3);

        dst_c1 = src_blend * (src_c1 - dst_c1) + dst_c1;
        dst_c2 = src_blend * (src_c2 - dst_c2) + dst_c2;
        dst_c3 = src_blend * (src_c3 - dst_c3) + dst_c3;

    } else {
        if (!haveMask && !haveOpacity) {
            memcpy(dst, src, 4 * float_v::size());
            return;
        } else {
            // opacity has changed the alpha of the source,
            // so we can't just memcpy the bytes
            dst_c1 = src_c1;
            dst_c2 = src_c2;
            dst_c3 = src_c3;
        }
    }

    KoStreamedMath::write_channels_32(dst, new_alpha, dst_c1, dst_c2, dst_c3);
}

template<typename T>
ALWAYS_INLINE T calculateZeroFlowAlpha(T srcAlpha, T dstAlpha, T normCoeff) {
    return srcAlpha + dstAlpha - srcAlpha * dstAlpha * normCoeff;
}

template<bool haveMask, bool src_aligned>
ALWAYS_INLINE void compositeVectorAlphaDarken(const quint8 *src, quint8 *dst, const quint8 *mask, float opacity)
{
    float_v src_alpha;
    float_v dst_alpha;

    // we don't use directly passed value
    Q_UNUSED(opacity);

    const float hackOpacity = 1.0;
    const float hackAverageOpacity = 1.0;
    const float hackFlow = 1.0;

    // instead we use value calculated by ParamsWrapper
    opacity = hackOpacity;
    float_v opacity_vec(255.0 * opacity);

    float_v average_opacity_vec(255.0 * hackAverageOpacity);
    float_v flow_norm_vec(hackFlow);


    float_v uint8MaxRec2((float)1.0 / (255.0 * 255.0));
    float_v uint8MaxRec1((float)1.0 / 255.0);
    float_v uint8Max((float)255.0);
    float_v zeroValue(Vc::Zero);


    float_v msk_norm_alpha;
    src_alpha = KoStreamedMath::template fetch_alpha_32<src_aligned>(src);

    if (haveMask) {
        float_v mask_vec = KoStreamedMath::fetch_mask_8(mask);
        msk_norm_alpha = src_alpha * mask_vec * uint8MaxRec2;
    } else {
        msk_norm_alpha = src_alpha * uint8MaxRec1;
    }

    dst_alpha = KoStreamedMath::template fetch_alpha_32<true>(dst);
    src_alpha = msk_norm_alpha * opacity_vec;

    float_m empty_dst_pixels_mask = dst_alpha == zeroValue;

    float_v src_c1;
    float_v src_c2;
    float_v src_c3;

    float_v dst_c1;
    float_v dst_c2;
    float_v dst_c3;

    KoStreamedMath::template fetch_colors_32<src_aligned>(src, src_c1, src_c2, src_c3);

    bool srcAlphaIsZero = (src_alpha == zeroValue).isFull();
    if (srcAlphaIsZero) return;

    bool dstAlphaIsZero = empty_dst_pixels_mask.isFull();

    float_v dst_blend = src_alpha * uint8MaxRec1;

    bool srcAlphaIsUnit = (src_alpha == uint8Max).isFull();

    if (dstAlphaIsZero) {
        dst_c1 = src_c1;
        dst_c2 = src_c2;
        dst_c3 = src_c3;
    } else if (srcAlphaIsUnit) {
        bool dstAlphaIsUnit = (dst_alpha == uint8Max).isFull();
        if (dstAlphaIsUnit) {
            memcpy(dst, src, 4 * float_v::size());
            return;
        } else {
            dst_c1 = src_c1;
            dst_c2 = src_c2;
            dst_c3 = src_c3;
        }
    } else if (empty_dst_pixels_mask.isEmpty()) {
        KoStreamedMath::template fetch_colors_32<true>(dst, dst_c1, dst_c2, dst_c3);
        dst_c1 = dst_blend * (src_c1 - dst_c1) + dst_c1;
        dst_c2 = dst_blend * (src_c2 - dst_c2) + dst_c2;
        dst_c3 = dst_blend * (src_c3 - dst_c3) + dst_c3;
    } else {
        KoStreamedMath::template fetch_colors_32<true>(dst, dst_c1, dst_c2, dst_c3);
        dst_c1(empty_dst_pixels_mask) = src_c1;
        dst_c2(empty_dst_pixels_mask) = src_c2;
        dst_c3(empty_dst_pixels_mask) = src_c3;

        float_m not_empty_dst_pixels_mask = !empty_dst_pixels_mask;

        dst_c1(not_empty_dst_pixels_mask) = dst_blend * (src_c1 - dst_c1) + dst_c1;
        dst_c2(not_empty_dst_pixels_mask) = dst_blend * (src_c2 - dst_c2) + dst_c2;
        dst_c3(not_empty_dst_pixels_mask) = dst_blend * (src_c3 - dst_c3) + dst_c3;
    }

    float_v fullFlowAlpha;

    if (hackAverageOpacity > opacity) {
        float_m fullFlowAlpha_mask = average_opacity_vec > dst_alpha;

        if (fullFlowAlpha_mask.isEmpty()) {
            fullFlowAlpha = dst_alpha;
        } else {
            float_v reverse_blend = dst_alpha / average_opacity_vec;
            float_v opt1 = (average_opacity_vec - src_alpha) * reverse_blend + src_alpha;
            fullFlowAlpha(!fullFlowAlpha_mask) = dst_alpha;
            fullFlowAlpha(fullFlowAlpha_mask) = opt1;
        }
    } else {
        float_m fullFlowAlpha_mask = opacity_vec > dst_alpha;

        if (fullFlowAlpha_mask.isEmpty()) {
            fullFlowAlpha = dst_alpha;
        } else {
            float_v opt1 = (opacity_vec - dst_alpha) * msk_norm_alpha + dst_alpha;
            fullFlowAlpha(!fullFlowAlpha_mask) = dst_alpha;
            fullFlowAlpha(fullFlowAlpha_mask) = opt1;
        }
    }

    if (hackFlow == 1.0) {
        dst_alpha = fullFlowAlpha;
    } else {
        float_v zeroFlowAlpha = calculateZeroFlowAlpha(src_alpha, dst_alpha, uint8MaxRec1);
        dst_alpha = (fullFlowAlpha - zeroFlowAlpha) * flow_norm_vec + zeroFlowAlpha;
    }

    KoStreamedMath::write_channels_32(dst, dst_alpha, dst_c1, dst_c2, dst_c3);
}

template<typename T>
ALWAYS_INLINE
T pow2(const T& x) {
    return x * x;
}


const int brushSize = 1000;
const qreal fade = 0.5;
const qreal angle = 0.18;
const qreal softness = 1.0;

const qreal xcoef = 2.0 / brushSize;
const qreal ycoef = 2.0 / brushSize;
const qreal xfadecoef = 2.0 / (fade * brushSize);
const qreal yfadecoef = 2.0 / (fade * brushSize);
const qreal transformedFadeX = softness * xfadecoef;
const qreal transformedFadeY = softness * yfadecoef;


template <bool useSmoothing, bool noFading>
ALWAYS_INLINE void processBrushMaskLine(float* buffer, int width, float y, float cosa, float sina,
                                        float centerX, float centerY)
{
    float y_ = y - centerY;
    float sinay_ = sina * y_;
    float cosay_ = cosa * y_;

    float* bufferPointer = buffer;

    float_v currentIndices = float_v::IndexesFromZero();

    float_v increment((float)float_v::size());
    float_v vCenterX(centerX);

    float_v vCosa(cosa);
    float_v vSina(sina);
    float_v vCosaY_(cosay_);
    float_v vSinaY_(sinay_);

    float_v vXCoeff(xcoef);
    float_v vYCoeff(ycoef);

    float_v vTransformedFadeX(transformedFadeX);
    float_v vTransformedFadeY(transformedFadeY);

    float_v vOne(Vc::One);

    for (int i=0; i < width; i+= float_v::size()){

        float_v x_ = currentIndices - vCenterX;

        float_v xr = x_ * vCosa - vSinaY_;
        float_v yr = x_ * vSina + vCosaY_;

        float_v n = pow2(xr * vXCoeff) + pow2(yr * vYCoeff);
        float_m outsideMask = n > vOne;

        if (!outsideMask.isFull()) {

            if (noFading) {
                float_v vFade(Vc::Zero);
                vFade(outsideMask) = vOne;
                vFade.store(bufferPointer, Vc::Aligned);
            } else {
                if (useSmoothing) {
                    xr = Vc::abs(xr) + vOne;
                    yr = Vc::abs(yr) + vOne;
                }

                float_v vNormFade = pow2(xr * vTransformedFadeX) + pow2(yr * vTransformedFadeY);

                //255 * n * (normeFade - 1) / (normeFade - n)
                float_v vFade = n * (vNormFade - vOne) / (vNormFade - n);

                // Mask in the inner circle of the mask
                float_m mask = vNormFade < vOne;
                vFade.setZero(mask);

                // Mask out the outer circle of the mask
                vFade(outsideMask) = vOne;

                vFade.store(bufferPointer, Vc::Aligned);
            }
        } else {
            // Mask out everything outside the circle
            vOne.store(bufferPointer, Vc::Aligned);
        }

        currentIndices = currentIndices + increment;

        bufferPointer += float_v::size();
    }
}

void testBrushMaskSpeed()
{
    const int numPixels = brushSize * brushSize;
    const int pixelAlignment = 32;
    const int srcAlignmentShift = 0;
    const int pixelSize = sizeof(float);

    int error = 0;

    void *srcPtr = 0;
    error = MEMALIGN_ALLOC(&srcPtr, pixelAlignment, numPixels * pixelSize + srcAlignmentShift);
    if (error) {
        qDebug() << "1";
        qFatal("posix_memalign failed: %d", error);
    }

    QElapsedTimer t;
    t.start();


    const float sina = std::sin(angle);
    const float cosa = std::cos(angle);

    for (int i = 0; i < 1000; i++) {

        float *src = reinterpret_cast<float*>(srcPtr);
        for (int y = 0; y < brushSize; y++) {
            processBrushMaskLine<true, false>(src, brushSize, float(y), cosa, sina,
                                              500.0, 500.0);

            src += brushSize;
        }
    }

    qDebug() << "Brush test time:" << t.elapsed() << "ms";

    MEMALIGN_FREE(srcPtr);
}

struct AlphaDarkenOp {
    void operator() (const quint8 *src, quint8 *dst, const quint8 *mask, float opacity) {
        compositeVectorAlphaDarken<true, true>(src, dst, mask, opacity);
    }

    static constexpr const char* name = "Alpha Darken";
};

struct OverOp {
    void operator() (const quint8 *src, quint8 *dst, const quint8 *mask, float opacity) {
        compositeVectorOver<true, true>(src, dst, mask, opacity);
    }

    static constexpr const char* name = "Over";
};


template <typename Func>
void testCompositionSpeed()
{
    const int numPixels = 100024000;
    const int pixelAlignment = 32;
    const int srcAlignmentShift = 0;
    const int dstAlignmentShift = 0;
    const int maskAlignmentShift = 0;
    const int pixelSize = 4;

    int error = 0;

    void *srcPtr = 0;
    error = MEMALIGN_ALLOC(&srcPtr, pixelAlignment, numPixels * pixelSize + srcAlignmentShift);
    if (error) {
        qFatal("posix_memalign failed: %d", error);
    }

    void *dstPtr = 0;
    error = MEMALIGN_ALLOC(&dstPtr, pixelAlignment, numPixels * pixelSize + dstAlignmentShift);
    if (error) {
        qFatal("posix_memalign failed: %d", error);
    }

    void *maskPtr = 0;
    error = MEMALIGN_ALLOC(&maskPtr, pixelAlignment, numPixels + maskAlignmentShift);
    if (error) {
        qFatal("posix_memalign failed: %d", error);
    }

    const int pixelsPerBlock = float_v::Size;
    const int numBlocks = numPixels / pixelsPerBlock;

    QElapsedTimer t;
    t.start();

    for (int j = 0; j < 50; j++) {

        quint8 *src = reinterpret_cast<quint8*>(srcPtr);
        quint8 *dst = reinterpret_cast<quint8*>(dstPtr);
        quint8 *mask = reinterpret_cast<quint8*>(maskPtr);

        for (int i = 0; i < numBlocks; i++) {
            Func()(src, dst, mask, 1.0);
            src += pixelsPerBlock * pixelSize;
            dst += pixelsPerBlock * pixelSize;
            mask += pixelsPerBlock;
        }
    }

    qDebug() << "Composition test time:" << Func::name << t.elapsed() << "ms";

    MEMALIGN_FREE(srcPtr);
    MEMALIGN_FREE(dstPtr);
    MEMALIGN_FREE(maskPtr);
}

int main(int argc, char *argv[])
{
    Q_UNUSED(argc);
    Q_UNUSED(argv);

    testCompositionSpeed<OverOp>();
    testCompositionSpeed<AlphaDarkenOp>();
    testBrushMaskSpeed();
    return 0;
}
