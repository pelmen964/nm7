using System;
using System.IO;
using System.Text;
using AngouriMath;
using AngouriMath.Extensions;

namespace nm7
{
    public class Task
    {
        /*
         *  левосторонних и правосторонних прямоугольников, трапеций и Симпсона и,
         *  по желанию, один из двух дополнительных – Чебышева
         *  или Гаусса.
         *
         */

        private readonly string _outFileName = "TaskOutput.txt";
        private readonly Entity _func;
        private readonly int _numMethod;
        private readonly int _gridType;
        private readonly int _numberOfintervals;
        private readonly double[] _borders;
        private readonly double[] _xGridNodes;
        private readonly double[] _yGridNodes;
        private readonly double _eps;
        private double _integralSum;
        private int _iterationCount;
        private double _pression;
        private double[] _abscissasOfIntegrationPoints;
        private double[] _gaussCoeffs;

        public Task(string inputFileName)
        {
            using (StreamReader fstream = new StreamReader(inputFileName))
            {
                _numMethod = int.Parse(fstream.ReadLine() ?? "0");
                _gridType = int.Parse(fstream.ReadLine() ?? "0");
                _numberOfintervals = int.Parse(fstream.ReadLine() ?? "0");

                Converter<string, double> parseStringToDouble = str => double.Parse(str);
                _borders = Array.ConvertAll(fstream.ReadLine().Split(), parseStringToDouble);
                _xGridNodes = Array.ConvertAll(fstream.ReadLine().Split(), parseStringToDouble);
                _func = fstream.ReadLine();
                _yGridNodes = Array.ConvertAll(fstream.ReadLine().Split(), parseStringToDouble);
                _eps = double.Parse(fstream.ReadLine() ?? "0");
            }
        }

        public double F(double x, uint der = 0)
        {
            var localExpr = _func;
            var localDer = der;
            while (localDer > 0)
            {
                localExpr = localExpr.Differentiate("x");
                localDer -= 1;
            }

            return (double)(localExpr.Substitute("x", x).EvalNumerical());
        }

        public double IntegrateWithLefthanded()
        {
            _integralSum = 0;
            _pression = 0;
            switch (_gridType)
            {
                case 1: // сетка равномерная
                    double lenOfInterval = (_borders[1] - _borders[0]) / _numberOfintervals;
                    for (double i = _borders[0]; i < _borders[1]; i += lenOfInterval)
                    {
                        _integralSum += F(i) * lenOfInterval;
                        _pression += F(i+lenOfInterval/2, 1) * lenOfInterval * lenOfInterval / 2;
                    }

                    return _integralSum;

                case 2:
                    for (int i = 1; i < _xGridNodes.Length; i++)
                    {
                        _integralSum += _yGridNodes[i-1] * (_xGridNodes[i] - _xGridNodes[i - 1]);
                        _pression += F(_yGridNodes[i-1], 1) * (_xGridNodes[i] - _xGridNodes[i - 1]) * (_xGridNodes[i] - _xGridNodes[i - 1]) / 2;
                    }

                    return _integralSum;
                case 3: // Динамическая сетка (просто жеееесть)
                    double internalStep = _eps;
                    double currentX = _borders[0] + internalStep;
                    _integralSum += F(_borders[0]) * internalStep;
                    while (currentX < _borders[1])
                    {
                        _integralSum += F(currentX) * internalStep;
                        // Расчитываю по формуле 9.33 из большой методички стр. 169
                        internalStep = (_borders[1] - _borders[0]) / 2 * internalStep * F(currentX, 1);
                        _pression += (_borders[1] - _borders[0]) / 2 * internalStep * F(currentX, 1);
                    }

                    _integralSum += F(_borders[1]) * internalStep;

                    return _integralSum;
                default:
                    throw new ArgumentException();
            }
        }
        public double IntegrateWithRighthanded()
        {
            _integralSum = 0;
            _pression = 0;
            switch (_gridType)
            {
                case 1: // сетка равномерная
                    double lenOfInterval = (_borders[1] - _borders[0]) / _numberOfintervals;
                    for (double i = _borders[0]+lenOfInterval; i <= _borders[1]; i += lenOfInterval)
                    {
                        _integralSum += F(i) * lenOfInterval;
                        _pression += F(i+lenOfInterval/2, 1) * lenOfInterval * lenOfInterval / 2;
                    }

                    return _integralSum;

                case 2:
                    for (int i = 1; i < _xGridNodes.Length; i++)
                    {
                        _integralSum += _yGridNodes[i] * (_xGridNodes[i] - _xGridNodes[i - 1]);
                        _pression += F(_yGridNodes[i-1], 1) * (_xGridNodes[i] - _xGridNodes[i - 1]) * (_xGridNodes[i] - _xGridNodes[i - 1]) / 2
                    }

                    return _integralSum;
                case 3:
                    double internalStep = _eps;
                    double currentX = _borders[0] + internalStep;
                    while (currentX <= _borders[1])
                    {
                        _integralSum += F(currentX) * internalStep;
                        // Расчитываю по формуле 9.33 из большой методички стр. 169
                        internalStep = (_borders[1] - _borders[0]) / 2 * internalStep * F(currentX, 1);
                        _pression = (_borders[1] - _borders[0]) / 2 * internalStep * F(currentX, 1);
                    }

                    _integralSum += F(_borders[1]) * internalStep;

                    return _integralSum;
                default:
                    throw new ArgumentException();
            }
        }
        public double 
        
    }
}