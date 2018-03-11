using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace Noise
{
    /// <summary>
    /// Implementation of simplex noise with analytic derivatives.
    /// Ported from Brian Sharpes GLSL noise library:
    /// https://github.com/BrianSharpe/GPU-Noise-Lib/blob/master/gpu_noise_lib.glsl
    /// </summary>
    public class SimplexNoise : NoiseModule
    {
        /// <summary>
        /// 2D version of simplex noise.
        /// </summary>
        /// <param name="position"></param>
        /// <param name="derivative"></param>
        /// <returns></returns>
        public override float Evaluate(Vector2 position, out Vector2 derivative)
        {
            throw new NotImplementedException();
        }


#region Vector Swizzling
        // VECTOR 3s
        static Vector3 XZY(Vector3 vec)
        {
            return new Vector3(vec.x, vec.z, vec.y);
        }

        static Vector3 YXZ(Vector3 vec)
        {
            return new Vector3(vec.y, vec.x, vec.z);
        }

        static Vector3 YZX(Vector3 vec)
        {
            return new Vector3(vec.y, vec.z, vec.x);
        }

        static Vector3 ZXY(Vector3 vec)
        {
            return new Vector3(vec.z, vec.x, vec.y);
        }

        static Vector3 ZYX(Vector3 vec)
        {
            return new Vector3(vec.z, vec.y, vec.x);
        }


#endregion

#region GLSL Utility Functions
        
        Vector4 InverseSqrt(Vector4 a)
        {
            return new Vector4(
                1f / Mathf.Sqrt(a.x),
                1f / Mathf.Sqrt(a.y),
                1f / Mathf.Sqrt(a.z),
                1f / Mathf.Sqrt(a.w));
        }

        Vector3 InverseSqrt(Vector3 a)
        {
            return new Vector4(
                1f / Mathf.Sqrt(a.x),
                1f / Mathf.Sqrt(a.y),
                1f / Mathf.Sqrt(a.z));
        }

        static Vector4 Fract(Vector4 a)
        {
            return new Vector4(
                a.x - Mathf.Floor(a.x),
                a.y - Mathf.Floor(a.y),
                a.z - Mathf.Floor(a.z),
                a.w - Mathf.Floor(a.w));
        }

        static Vector4 Multi(Vector4 a, Vector4 b)
        {
            return new Vector4(a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w);
        }

        static Vector3 Divide(Vector3 a, Vector3 b)
        {
            return new Vector3(a.x / b.x, a.y / b.y, a.z / b.z);
        }

        static Vector3 Multi(Vector3 a, Vector3 b)
        {
            return new Vector3(a.x * b.x, a.y * b.y, a.z * b.z);
        }

        static Vector2 Min(Vector2 v1, Vector2 v2)
        {
            return new Vector2(
                v1.x < v2.x ? v1.x : v2.x,
                v1.y < v2.y ? v1.y : v2.y);
        }

        static Vector3 Min(Vector3 v1, Vector3 v2)
        {
            return new Vector3(
                v1.x < v2.x ? v1.x : v2.x,
                v1.y < v2.y ? v1.y : v2.y,
                v1.z < v2.z ? v1.z : v2.z);
        }

        static Vector4 Min(Vector4 v1, Vector4 v2)
        {
            return new Vector4(
                v1.x < v2.x ? v1.x : v2.x,
                v1.y < v2.y ? v1.y : v2.y,
                v1.z < v2.z ? v1.z : v2.z,
                v1.w < v2.w ? v1.w : v2.w);
        }

        static Vector2 Max(Vector2 v1, Vector2 v2)
        {
            return new Vector2(
                v1.x > v2.x ? v1.x : v2.x,
                v1.y > v2.y ? v1.y : v2.y);
        }

        static Vector3 Max(Vector3 v1, Vector3 v2)
        {
            return new Vector3(
                v1.x > v2.x ? v1.x : v2.x,
                v1.y > v2.y ? v1.y : v2.y,
                v1.z > v2.z ? v1.z : v2.z);
        }

        static Vector4 Max(Vector4 v1, Vector4 v2)
        {
            return new Vector4(
                v1.x > v2.x ? v1.x : v2.x,
                v1.y > v2.y ? v1.y : v2.y,
                v1.z > v2.z ? v1.z : v2.z,
                v1.w > v2.w ? v1.w : v2.w);
        }

        static Vector2 Step(Vector2 x, Vector2 edge)
        {
            return new Vector2(
                x.x < edge.x ? 0f : 1f,
                x.y < edge.y ? 0f : 1f);
        }

        static Vector3 Step(Vector3 x, Vector3 edge)
        {
            return new Vector3(
                x.x < edge.x ? 0f : 1f,
                x.y < edge.y ? 0f : 1f,
                x.z < edge.z ? 0f : 1f);
        }

        static Vector4 Step(Vector4 x, Vector4 edge)
        {
            return new Vector4(
                x.x < edge.x ? 0f : 1f,
                x.y < edge.y ? 0f : 1f,
                x.z < edge.z ? 0f : 1f,
                x.w < edge.w ? 0f : 1f);
        }

        static float Dot(Vector3 a, Vector3 b)
        {
            return Vector3.Dot(a, b);
        }

        static Vector4 Vec4(float f)
        {
            return new Vector4(f, f, f, f);
        }

        static Vector4 Vec4(float x, float y, float z, float w)
        {
            return new Vector4(x, y, z, w);
        }

        static Vector3 Vec3(float f)
        {
            return new Vector3(f, f, f);
        }

        static Vector3 Vec3(float x, float y, float z)
        {
            return new Vector3(x, y, z);
        }

        static Vector2 Vec2(float x, float y)
        {
            return new Vector2(x, y);
        }

        static Vector4 Floor(Vector4 vec)
        {
            return new Vector4(Mathf.Floor(vec.x), Mathf.Floor(vec.y), Mathf.Floor(vec.z), Mathf.Floor(vec.w));
        }

        static Vector3 Floor(Vector3 vec)
        {
            return new Vector3(Mathf.Floor(vec.x), Mathf.Floor(vec.y), Mathf.Floor(vec.z));
        }

        static Vector2 Floor(Vector2 vec)
        {
            return new Vector2(Mathf.Floor(vec.x), Mathf.Floor(vec.y));
        }

        static Vector3 Mix(Vector3 a, Vector3 b, float t)
        {
            return Vector3.Lerp(a, b, t);
        }

        static Vector4 Mix(Vector4 a, Vector4 b, Vector4 t)
        {
            return new Vector4(
                Mathf.Lerp(a.x, b.x, t.x),
                Mathf.Lerp(a.y, b.y, t.y),
                Mathf.Lerp(a.z, b.z, t.z),
                Mathf.Lerp(a.w, b.w, t.w));
        }

        #endregion

        static void GetCornerVectors(Vector3 p, out Vector3 pi, out Vector3 pi1, out Vector3 pi2, out Vector4 vx, out Vector4 vy, out Vector4 vz)
        {
            //
            //	Simplex math from Stefan Gustavson's and Ian McEwan's work at...
            //	http://github.com/ashima/webgl-noise
            //
            //	simplex math constants
            const float SKEWFACTOR = 1.0f / 3.0f;
            const float UNSKEWFACTOR = 1.0f / 6.0f;
            const float SIMPLEX_CORNER_POS = 0.5f;
            const float SIMPLEX_PYRAMID_HEIGHT = 0.70710678118654752440084436210485f;    // sqrt( 0.5 )	height of simplex pyramid.

            p *= SIMPLEX_PYRAMID_HEIGHT;        // scale space so we can have an approx feature size of 1.0  ( optional )

            //	Find the vectors to the corners of our simplex pyramid
            pi = Floor(p + Vec3(Dot(p, Vec3(SKEWFACTOR))));
            Vector3 x0 = p - pi + Vec3(Dot(pi, Vec3(UNSKEWFACTOR)));
            Vector3 g = Step(YZX(x0), x0);
            Vector3 l = Vec3(1.0f) - g;
            pi1 = Min(g, ZXY(l));
            pi2 = Max(g, ZXY(l));
            Vector3 x1 = x0 - pi1 + Vec3(UNSKEWFACTOR);
            Vector3 x2 = x0 - pi2 + Vec3(SKEWFACTOR);
            Vector3 x3 = x0 - Vec3(SIMPLEX_CORNER_POS);

            //	pack them into a parallel-friendly arrangement
            vx = Vec4(x0.x, x1.x, x2.x, x3.x);
            vy = Vec4(x0.y, x1.y, x2.y, x3.y);
            vz = Vec4(x0.z, x1.z, x2.z, x3.z);
        }

        static void Hash32(Vector3 gridcell, Vector3 v1Mask, Vector3 v2Mask, out Vector4 hash0, out Vector4 hash1, out Vector4 hash2)
        {
            //    gridcell is assumed to be an integer coordinate

            //	TODO: 	these constants need tweaked to find the best possible noise.
            //			probably requires some kind of brute force computational searching or something...
            Vector2 OFFSET = Vec2(50.0f, 161.0f);
            const float DOMAIN = 69.0f;
            Vector3 SOMELARGEFLOATS = Vec3(635.298681f, 682.357502f, 668.926525f);
            Vector3 ZINC = Vec3(48.500388f, 65.294118f, 63.934599f);

            //	truncate the domain
            gridcell = gridcell - Floor(gridcell * (1.0f / DOMAIN)) * DOMAIN;
            Vector3 gridcell_inc1 = Step(gridcell, Multi(Vec3(DOMAIN - 1.5f), gridcell + Vec3(1.0f)));

            //	compute x*x*y*y for the 4 corners
            Vector4 p = Vec4(gridcell.x, gridcell.y, gridcell_inc1.x, gridcell_inc1.y) + Vec4(OFFSET.x, OFFSET.y, OFFSET.x, OFFSET.y);
            p.x *= p.x;
            p.y *= p.y;
            p.z *= p.z;
            p.w *= p.w;
            Vector4 V1xy_V2xy = Mix(Vec4(p.x, p.y, p.x, p.y), Vec4(p.z, p.w, p.z, p.w), Vec4(v1Mask.x, v1Mask.y, v2Mask.x, v2Mask.y));     //	apply mask for v1 and v2
            p = Multi(Vec4(p.x, V1xy_V2xy.x, V1xy_V2xy.z, p.z), Vec4(p.y, V1xy_V2xy.y, V1xy_V2xy.w, p.w));

            //	get the lowz and highz mods
            Vector3 lowz_mods = Divide(Vec3(1.0f), (SOMELARGEFLOATS + Multi(Vec3(gridcell.z), ZINC)));
            Vector3 highz_mods = Divide(Vec3(1.0f), (SOMELARGEFLOATS + Multi(Vec3(gridcell_inc1.z), ZINC)));

            //	apply mask for v1 and v2 mod values
            v1Mask = (v1Mask.z < 0.5) ? lowz_mods : highz_mods;
            v2Mask = (v2Mask.z < 0.5) ? lowz_mods : highz_mods;

            //	compute the final hash
            hash0 = Fract(Multi(p, Vec4(lowz_mods.x, v1Mask.x, v2Mask.x, highz_mods.x)));
            hash1 = Fract(Multi(p, Vec4(lowz_mods.y, v1Mask.y, v2Mask.y, highz_mods.y)));
            hash2 = Fract(Multi(p, Vec4(lowz_mods.z, v1Mask.z, v2Mask.z, highz_mods.z)));
        }

        /// <summary>
        /// 3D version of simplex noise.
        /// </summary>
        /// <param name="position"></param>
        /// <param name="derivative"></param>
        /// <returns></returns>
        public override float Evaluate(Vector3 position, out Vector3 derivative)
        {
            //	calculate the simplex vector and index math
            Vector3 pi;
            Vector3 pi1;
            Vector3 pi2;
            Vector4 vx;
            Vector4 vy;
            Vector4 vz;
            GetCornerVectors(position, out pi, out pi1, out pi2, out vx, out vy, out vz);

            //	generate the random vectors
            //	( various hashing methods listed in order of speed )
            Vector4 hash_0;
            Vector4 hash_1;
            Vector4 hash_2;
            Hash32(pi, pi1, pi2, out hash_0, out hash_1, out hash_2);
            //SGPP_hash_3D( Pi, Pi_1, Pi_2, hash_0, hash_1, hash_2 );
            hash_0 -= Vec4(0.49999f);
            hash_1 -= Vec4(0.49999f);
            hash_2 -= Vec4(0.49999f);

            //	normalize random gradient vectors
            Vector4 norm = InverseSqrt(Multi(hash_0, hash_0) + Multi(hash_1, hash_1) + Multi(hash_2, hash_2));
            hash_0 = Multi(hash_0, norm);
            hash_1 = Multi(hash_1, norm);
            hash_2 = Multi(hash_2, norm);

            //	evaluate gradients
            Vector4 grad_results = Multi(hash_0, vx) + Multi(hash_1, vy) + Multi(hash_2, vz);

            //	evaluate the surflet f(x)=(0.5-x*x)^3
            Vector4 m = Multi(vx, vx) + Multi(vy, vy) + Multi(vz, vz);
            m = Max(Vec4(0.5f) - m, Vec4(0f));      //	The 0.5 here is SIMPLEX_PYRAMID_HEIGHT^2
            Vector4 m2 = Multi(m, m);
            Vector4 m3 = Multi(m, m2);

            //	calc the deriv
            Vector4 temp = Multi(-6.0f * m2, grad_results);
            float xderiv = Dot(temp, vx) + Dot(m3, hash_0);
            float yderiv = Dot(temp, vy) + Dot(m3, hash_1);
            float zderiv = Dot(temp, vz) + Dot(m3, hash_2);

            const float FINAL_NORMALIZATION = 37.837227241611314102871574478976f;    //	scales the final result to a strict 1.0->-1.0 range
            derivative = new Vector3(xderiv, yderiv, zderiv);
            //	sum with the surflet and return
            return Dot(m3, grad_results) * FINAL_NORMALIZATION;
        }
    }
}