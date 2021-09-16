using System;

namespace BPSE
{
    class VonMises
    {
        /// <summary>
        /// フォン・ミーゼス分布を初期化します。
        /// </summary>
        /// <param name="m">集中度パラメター</param>
        /// <param name="iterate">0次第1種変形ベッセル関数の近似計算の反復回数</param>
        public VonMises(double theta0, double m, int iterate = 100)
        {
            _theta0 = theta0;
            _m = m;
            _I0 = Bessel_I0(m, iterate);
        }
        private readonly double _theta0;
        private readonly double _m;
        private readonly double _I0;
        /// <summary>
        /// フォン・ミーゼス分布の確率密度を取得します。　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　
        /// </summary>
        /// <param name="theta">微小範囲[rad]</param>
        /// <returns></returns>
        public double Get(double theta)
        {
            //return 1.0 / (2.0 * Math.PI * _I0) * Math.Exp(_m * Math.Cos(theta - _theta0));
            return Math.Exp(_m * Math.Cos(theta - _theta0)) / (2.0 * Math.PI * _I0);
        }
        /// <summary>
        /// フォン・ミーゼス分布の確率密度を取得します。　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　
        /// </summary>
        /// <param name="theta">微小範囲[rad]</param>
        /// <param name="theta0">分布の平均に相当 [rad]</param>
        /// <param name="m">集中度パラメター</param>
        /// <param name="iterate">0次第1種変形ベッセル関数の近似計算の反復回数</param>
        /// <returns></returns>
        public static double Get(double theta, double theta0, double m, int iterate = 100)
        {
            double I0 = Bessel_I0(m, iterate);
            //return 1.0 / (2.0 * Math.PI * I0) * Math.Exp(m * Math.Cos(theta - theta0));
            return Math.Exp(m * Math.Cos(theta - theta0)) / (2.0 * Math.PI * I0);
        }

        /// <summary>
        /// 0次第1種変形ベッセル関数の値を取得します。
        /// </summary>
        /// <param name="x">x の値。</param>
        /// <param name="iterate">近似計算の反復回数。</param>
        /// <returns>計算結果。</returns>
        public static double Bessel_I0(double x, int iterate = 100)
        {
            if (x < 0.0)
                throw new ArgumentOutOfRangeException(nameof(x));

            if (iterate < 0)
                throw new ArgumentOutOfRangeException(nameof(iterate));

            double y = 0.0;
            for (var i = 0; i < iterate; i++)
            {
#if true
                double y_ = 1.0;
                for (int n = 1; n <= i; n++)
                    y_ *= x / n / 2.0;
                if (y_ == 0.0) break;
                y += y_ * y_;
#else
                double f = InvertedFactorial(i);
                if (f == 0) break;
                double y2 = Math.Pow(x / 2.0, 2.0 * i);
                if (double.IsPositiveInfinity(y2)) break;
                y += f * f * y2;
#endif
            }

            return y;
        }

        /// <summary>
        /// n の階乗を逆数で計算します。
        /// </summary>
        /// <param name="n">階乗の整数パラメータ。</param>
        /// <returns>n の階乗の逆数の値。</returns>
        private static double InvertedFactorial(int n)
        {
            if (n < 0)
                throw new ArgumentOutOfRangeException(nameof(n));

            var y = 1.0;

            while (n > 1)
            {
                y *= 1.0 / n;
                n--;
            }

            return y;
        }
    }
}
