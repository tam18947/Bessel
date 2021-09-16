using System;

namespace Bessel
{
    class Program
    {
        static void Main(string[] args)
        {
            Console.WriteLine("modified Bessel function of the first kind");

            for (int m = 0; m < 1000; m++)
            {
                double val = BPSE.VonMises.Bessel_I0(m * 0.1, 1000);
                Console.WriteLine(val);
            }

            Console.WriteLine("von Mises");
            BPSE.VonMises vonMises = new BPSE.VonMises(Math.PI / 4.0, 5);
            for (int i = 0; i < 360; i++)
            {
                double theta = i * Math.PI / 180;
                double val = BPSE.VonMises.Get(theta, Math.PI / 4.0, 5);
                double val2 = vonMises.Get(theta);
                Console.WriteLine(val + "\t" + (val - val2));
            }
        }
    }
}
