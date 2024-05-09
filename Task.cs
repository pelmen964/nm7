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
        private double _accuracy;
        private double[] _abscissasOfIntegrationPoints;
        private double[] _gaussCoeffs;
        private readonly int _funcType;

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
                _funcType = int.Parse(fstream.ReadLine() ?? "0");
                _yGridNodes = Array.ConvertAll(fstream.ReadLine().Split(), parseStringToDouble);
                _func = fstream.ReadLine();
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

        public string IntegrateWithLefthanded()
        {
            var sb = new StringBuilder();
            sb.AppendLine("Левые прямоуголники:");
            _integralSum = 0;
            _accuracy = 0;
            switch (_gridType)
            {
                case 1: // сетка равномерная
                    double lenOfInterval = (_borders[1] - _borders[0]) / _numberOfintervals;
                    for (double i = _borders[0]; i < _borders[1]; i += lenOfInterval)
                    {
                        _integralSum += F(i) * lenOfInterval;
                        _accuracy += F(i + lenOfInterval / 2, 1) * lenOfInterval * lenOfInterval / 2;
                    }

                    break;

                case 2: // неравномерная сетка
                    if (_funcType == 1) // Табличные значения
                    {
                        for (int i = 1; i < _xGridNodes.Length; i++)
                        {
                            _integralSum += _yGridNodes[i - 1] * (_xGridNodes[i] - _xGridNodes[i - 1]);
                            _accuracy += F(_xGridNodes[i - 1], 1) * (_xGridNodes[i] - _xGridNodes[i - 1]) *
                                (_xGridNodes[i] - _xGridNodes[i - 1]) / 2;
                        }
                    }
                    else // Аналитическая функция
                    {
                        for (int i = 1; i < _xGridNodes.Length; i++)
                        {
                            _integralSum += (F(_xGridNodes[i - 1]) * (_xGridNodes[i] - _xGridNodes[i - 1]));
                            _accuracy += F(_xGridNodes[i - 1], 1) * (_xGridNodes[i] - _xGridNodes[i - 1]) *
                                (_xGridNodes[i] - _xGridNodes[i - 1]) / 2;
                        }
                    }

                    break;

                case 3: // Динамическая сетка (просто жеееесть)
                    double internalStep = 0;
                    double currentX = _borders[0];
                    int cnt = 0;
                    while (currentX < _borders[1])
                    {
                        cnt++;
                        _integralSum += F(currentX) * internalStep;
                        internalStep = F(currentX, 1) * _eps;
                        if (internalStep == 0)
                        {
                            internalStep = _eps;
                        }

                        // Расчитываю по формуле 9.33 из большой методички стр. 169
                        _accuracy += F(currentX, 1) / 2 * internalStep * internalStep;
                        currentX += internalStep;
                    }

                    _integralSum += F(_borders[1]) * internalStep;
                    _accuracy += F(_borders[1], 1) / 2 * internalStep * internalStep;
                    sb.AppendLine("Iters: " + cnt);
                    break;

                default:
                    throw new ArgumentException();
            }

            sb.AppendLine("Sum: " + _integralSum);
            sb.AppendLine("Eps: " + _accuracy);
            return sb.ToString();
        }

        public string IntegrateWithRighthanded()
        {
            var sb = new StringBuilder();
            sb.AppendLine("Правые прямоуголники:");
            _integralSum = 0;
            _accuracy = 0;
            switch (_gridType)
            {
                case 1: // сетка равномерная
                    double lenOfInterval = (_borders[1] - _borders[0]) / _numberOfintervals;
                    for (double i = _borders[0] + lenOfInterval; i <= _borders[1]; i += lenOfInterval)
                    {
                        _integralSum += F(i) * lenOfInterval;
                        _accuracy += F(i, 1) * lenOfInterval * lenOfInterval / 2;
                    }

                    break;

                case 2: // неравномерная сетка
                    if (_funcType == 1) // Табличные значения
                    {
                        for (int i = 1; i < _xGridNodes.Length; i++)
                        {
                            _integralSum += _yGridNodes[i] * (_xGridNodes[i] - _xGridNodes[i - 1]);
                            _accuracy += F(_xGridNodes[i], 1) * (_xGridNodes[i] - _xGridNodes[i - 1]) *
                                (_xGridNodes[i] - _xGridNodes[i - 1]) / 2;
                        }
                    }
                    else // Аналитическая функция
                    {
                        for (int i = 1; i < _xGridNodes.Length; i++)
                        {
                            _integralSum += (F(_xGridNodes[i]) * (_xGridNodes[i] - _xGridNodes[i - 1]));
                            _accuracy += F(_xGridNodes[i], 1) * (_xGridNodes[i] - _xGridNodes[i - 1]) *
                                (_xGridNodes[i] - _xGridNodes[i - 1]) / 2;
                        }
                    }

                    break;

                case 3: // Динамическая сетка (просто жеееесть)
                    double internalStep = _eps;
                    double currentX = _borders[0] + internalStep;
                    _integralSum += F(_borders[0]) * internalStep;
                    int cnt = 0;
                    while (currentX < _borders[1])
                    {
                        cnt++;
                        _integralSum += F(currentX) * internalStep;
                        // Расчитываю по формуле 9.33 из большой методички стр. 169
                        _accuracy += F(currentX, 1) / 2 * internalStep * internalStep;
                        internalStep = F(currentX, 1) * _eps;
                        currentX += internalStep;
                    }

                    _integralSum += F(_borders[1]) * internalStep;
                    _accuracy += F(_borders[1], 1) / 2 * internalStep * internalStep;
                    sb.AppendLine("Iters: " + cnt);
                    break;

                default:
                    throw new ArgumentException();
            }

            sb.AppendLine("Sum: " + _integralSum);
            sb.AppendLine("Eps: " + _accuracy);
            return sb.ToString();
        }

        public string IntegrateWithTrapeze()
        {
            var sb = new StringBuilder();
            sb.AppendLine("Трапеции:");
            _integralSum = 0;
            _accuracy = 0;
            switch (_gridType)
            {
                case 1: // сетка равномерная
                    double lenOfInterval = (_borders[1] - _borders[0]) / _numberOfintervals;
                    for (double i = _borders[0]; i < _borders[1]; i += lenOfInterval)
                    {
                        _integralSum += (F(i) + F(i + lenOfInterval)) * lenOfInterval / 2;
                        // По формуле 9.12
                        _accuracy += Math.Pow(lenOfInterval, 3) / 12 * F(i + lenOfInterval / 2, 2);
                    }

                    break;

                case 2: // неравномерная сетка
                    if (_funcType == 1) // Табличные значения
                    {
                        for (int i = 1; i < _xGridNodes.Length; i++)
                        {
                            _integralSum += (_yGridNodes[i - 1] + _yGridNodes[i]) *
                                (_xGridNodes[i] - _xGridNodes[i - 1]) / 2;
                            _accuracy += Math.Pow(_xGridNodes[i] - _xGridNodes[i - 1], 3) / 12 *
                                         F(_xGridNodes[i] + _xGridNodes[i - 1] / 2, 2);
                        }
                    }
                    else // Аналитическая функция
                    {
                        for (int i = 1; i < _xGridNodes.Length; i++)
                        {
                            _integralSum += (F(_xGridNodes[i - 1]) + F(_xGridNodes[i])) *
                                (_xGridNodes[i] - _xGridNodes[i - 1]) / 2;
                            _accuracy += Math.Pow((_xGridNodes[i] - _xGridNodes[i - 1]), 3) / 12 *
                                         F(_xGridNodes[i] + _xGridNodes[i - 1] / 2, 2);
                        }
                    }

                    break;

                case 3: // Динамическая сетка
                    double internalStep = _eps;
                    double currentX = _borders[0] + internalStep;
                    int cnt = 0;
                    while (currentX < _borders[1])
                    {
                        cnt++;
                        _integralSum += (F(currentX) + F(currentX - internalStep)) * internalStep / 2;
                        _accuracy += Math.Pow(internalStep, 3) / 12 * F(currentX - internalStep / 2, 2);
                        internalStep = F(currentX, 1) * _eps;
                        currentX += internalStep;
                    }

                    _integralSum += (F(_borders[1]) + F(currentX - internalStep)) *
                        (_borders[1] - currentX + internalStep) / 2;
                    _accuracy += Math.Pow(currentX - _borders[1], 3) / 12 * F(currentX - internalStep / 2, 2);
                    sb.AppendLine("Iters: " + cnt);
                    break;

                default:
                    throw new ArgumentException();
            }

            sb.AppendLine("Sum: " + _integralSum);
            sb.AppendLine("Eps: " + _accuracy);
            return sb.ToString();
        }

        public string IntegrateWithSimpson()
        {
            var sb = new StringBuilder();
            sb.AppendLine("Симпсон:");
            _integralSum = 0;
            _accuracy = 0;
            double coeff, accuracyCoeff;
            switch (_gridType)
            {
                case 1: // Сетка равномерная
                    double lenOfInterval = (_borders[1] - _borders[0]) / _numberOfintervals;
                    double doubleLen = 2 * lenOfInterval;
                    coeff = 2 * lenOfInterval / 6;
                    accuracyCoeff = Math.Pow(lenOfInterval, 5);
                    for (double i = _borders[0] + doubleLen; i <= _borders[1]; i += doubleLen)
                    {
                        _integralSum += coeff * (F(i - doubleLen) + 4 * F(i - lenOfInterval) + F(i));
                        // По формуле 9.18
                        _accuracy += accuracyCoeff / 180 * F(i - lenOfInterval, 4);
                    }

                    break;

                case 3: // Динамическая сетка
                    double internalStep = _eps;
                    double currentX = _borders[0] + 2 * internalStep;
                    int cnt = 0;
                    while (currentX < _borders[1])
                    {
                        cnt++;
                        _integralSum += 2 * internalStep / 6 *
                                        (F(currentX - internalStep * 2) + 4 * F(currentX - internalStep) +
                                         F(currentX));

                        _accuracy += Math.Pow(internalStep, 5) / 180 * F(currentX - internalStep, 4);
                        internalStep = F(currentX, 1) * _eps;
                        currentX += 2 * internalStep;
                    }

                    double lastStep = 2 * internalStep;
                    double stepToBorder = _borders[1] - (currentX - lastStep);
                    _integralSum += 1 * stepToBorder / 6 *
                                    (F(currentX - lastStep) + 4 * F((_borders[1] + (currentX - lastStep)) / 2) +
                                     F(_borders[1]));
                    _accuracy += Math.Pow(lastStep / 2, 5) / 180 * F(_borders[1] - lastStep / 2, 4);
                    sb.AppendLine("Iters: " + cnt);
                    break;

                default:
                    throw new ArgumentException();
            }

            sb.AppendLine("Sum: " + _integralSum);
            sb.AppendLine("Eps: " + _accuracy);
            return sb.ToString();
        }

        public string IntegrateWithGauss()
        {
            var sb = new StringBuilder();
            sb.AppendLine("Гаус:");
            _integralSum = 0;
            _accuracy = 0;
            int n = 2; // Степень полинома
            double[] w = { 1.0, 1.0 }; // весовые коэффициенты из таблицы 9.1
            double[] x = { -0.577350269, +0.577350269 }; // значения аргумента из таблицы 9.1


            for (int i = 0; i < n; i++)
            {
                // По формуле 9.19
                _integralSum += 0.5 * (_borders[1] - _borders[0]) * w[i] *
                                F(0.5 * (_borders[0] + _borders[1]) + 0.5 * (_borders[1] - _borders[0]) * x[i]);
                sb.AppendLine($"abciss{i+1}: {x[i]}");
                sb.AppendLine($"A{i+1}: {w[i]}");
            }
            sb.AppendLine("Sum: " + _integralSum);
            _accuracy = Math.Pow((_borders[1] - _borders[0])/2, 5) / 135 * F((_borders[1] + _borders[0]), 4);
            sb.AppendLine("Eps: " + _accuracy);
            return sb.ToString();
        }
    }
}

