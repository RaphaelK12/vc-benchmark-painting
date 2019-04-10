/*
 *  Copyright (c) 2012 Dmitry Kazakov <dimula73@gmail.com>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; see the file COPYING.LIB.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 */

#ifndef __KOSTREAMED_MATH_H
#define __KOSTREAMED_MATH_H

#if defined _MSC_VER
// Lets shut up the "possible loss of data" and "forcing value to bool 'true' or 'false'
#pragma warning ( push )
#pragma warning ( disable : 4244 )
#pragma warning ( disable : 4800 )
#endif
#include <Vc/Vc>
#include <Vc/IO>
#if defined _MSC_VER
#pragma warning ( pop )
#endif

#include <stdint.h>
#include <KoAlwaysInline.h>
#include <iostream>

#define BLOCKDEBUG 0

#if !defined _MSC_VER
#pragma GCC diagnostic ignored "-Wcast-align"
#endif

#ifdef USE_SSE
using float_v = Vc::SSE::float_v;
#else
using float_v = Vc::float_v;
#endif
using float_m = typename float_v::mask_type;

struct KoStreamedMath {

using int_v = Vc::SimdArray<int, float_v::size()>;
using uint_v = Vc::SimdArray<unsigned int, float_v::size()>;


static inline quint8 round_float_to_uint(float value) {
    return quint8(value + float(0.5));
}

static inline quint8 lerp_mixed_u8_float(quint8 a, quint8 b, float alpha) {
    return round_float_to_uint(qint16(b - a) * alpha + a);
}

/**
 * Get a vector containing first float_v::size() values of mask.
 * Each source mask element is considered to be a 8-bit integer
 */
static inline float_v fetch_mask_8(const quint8 *data) {
    uint_v data_i(data);
    return Vc::simd_cast<float_v>(int_v(data_i));
}

/**
 * Get an alpha values from float_v::size() pixels 32-bit each
 * (4 channels, 8 bit per channel).  The alpha value is considered
 * to be stored in the most significant byte of the pixel
 *
 * \p aligned controls whether the \p data is fetched using aligned
 *            instruction or not.
 *            1) Fetching aligned data with unaligned instruction
 *               degrades performance.
 *            2) Fetching unaligned data with aligned instruction
 *               causes \#GP (General Protection Exception)
 */
template <bool aligned>
static inline float_v fetch_alpha_32(const quint8 *data) {
    uint_v data_i;
    if (aligned) {
        data_i.load((const quint32*)data, Vc::Aligned);
    } else {
        data_i.load((const quint32*)data, Vc::Unaligned);
    }

    return Vc::simd_cast<float_v>(int_v(data_i >> 24));
}

/**
 * Get color values from float_v::size() pixels 32-bit each
 * (4 channels, 8 bit per channel).  The color data is considered
 * to be stored in the 3 least significant bytes of the pixel.
 *
 * \p aligned controls whether the \p data is fetched using aligned
 *            instruction or not.
 *            1) Fetching aligned data with unaligned instruction
 *               degrades performance.
 *            2) Fetching unaligned data with aligned instruction
 *               causes \#GP (General Protection Exception)
 */
template <bool aligned>
static inline void fetch_colors_32(const quint8 *data,
                            float_v &c1,
                            float_v &c2,
                            float_v &c3) {
    int_v data_i;
    if (aligned) {
        data_i.load((const quint32*)data, Vc::Aligned);
    } else {
        data_i.load((const quint32*)data, Vc::Unaligned);
    }

    const quint32 lowByteMask = 0xFF;
    uint_v mask(lowByteMask);

    c1 = Vc::simd_cast<float_v>(int_v((data_i >> 16) & mask));
    c2 = Vc::simd_cast<float_v>(int_v((data_i >> 8)  & mask));
    c3 = Vc::simd_cast<float_v>(int_v( data_i        & mask));
}

/**
 * Pack color and alpha values to float_v::size() pixels 32-bit each
 * (4 channels, 8 bit per channel).  The color data is considered
 * to be stored in the 3 least significant bytes of the pixel, alpha -
 * in the most significant byte
 *
 * NOTE: \p data must be aligned pointer!
 */
static inline void write_channels_32(quint8 *data,
                                     float_v::AsArg alpha,
                                     float_v::AsArg c1,
                                     float_v::AsArg c2,
                                     float_v::AsArg c3) {
    /**
     * FIXME: make conversion float->int
     * use methematical rounding
     */

    const quint32 lowByteMask = 0xFF;

    // FIXME: Use single-instruction rounding + conversion
    //        The achieve that we need to implement Vc::iRound()

    uint_v mask(lowByteMask);
    uint_v v1 = uint_v(int_v(Vc::round(alpha))) << 24;
    uint_v v2 = (uint_v(int_v(Vc::round(c1))) & mask) << 16;
    uint_v v3 = (uint_v(int_v(Vc::round(c2))) & mask) <<  8;
    uint_v v4 = uint_v(int_v(Vc::round(c3))) & mask;
    v1 = v1 | v2;
    v3 = v3 | v4;
    (v1 | v3).store((quint32*)data, Vc::Aligned);
}

};

namespace KoStreamedMathFunctions {

template<int pixelSize>
ALWAYS_INLINE void clearPixel(quint8* dst);

template<>
ALWAYS_INLINE void clearPixel<4>(quint8* dst)
{
    quint32 *d = reinterpret_cast<quint32*>(dst);
    *d = 0;
}

template<>
ALWAYS_INLINE void clearPixel<16>(quint8* dst)
{
    quint64 *d = reinterpret_cast<quint64*>(dst);
    d[0] = 0;
    d[1] = 0;
}

template<int pixelSize>
ALWAYS_INLINE void copyPixel(const quint8 *src, quint8* dst);

template<>
ALWAYS_INLINE void copyPixel<4>(const quint8 *src, quint8* dst)
{
    const quint32 *s = reinterpret_cast<const quint32*>(src);
    quint32 *d = reinterpret_cast<quint32*>(dst);
    *d = *s;
}

template<>
ALWAYS_INLINE void copyPixel<16>(const quint8 *src, quint8* dst)
{
    const quint64 *s = reinterpret_cast<const quint64*>(src);
    quint64 *d = reinterpret_cast<quint64*>(dst);
    d[0] = s[0];
    d[1] = s[1];
}
}

#endif /* __KOSTREAMED_MATH_H */
