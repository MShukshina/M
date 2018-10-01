using System;
using System.IO;
using System.Threading.Tasks;


namespace Lab3_Чисмет
{
    class Program
    {
        /// <summary>
        /// Вывод информации об итерации
        /// </summary>
        public static void WriteInformation(int k, double[,] A, double[] PreviousValues, double[] b, int n, double delta)
        {
            Console.WriteLine("Номер итерации = " + k);
            Console.WriteLine("Норма невязки = " + NormVector(GradientR(A, PreviousValues, b, n), n));
            Console.WriteLine("Оценка нормы матрицы перехода = " + delta);
            Console.WriteLine("Оценка погрешности приближенного решения = " );
        }

        /// <summary>
        /// Вывод матрицы
        /// </summary>
        /// <param name="matrix"></param>
        /// <param name="n"></param>
        /// <param name="rounded"></param>
        public static void WriteMatrix(double[,] matrix, int n, bool rounded = true)
        {
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                    if (rounded)
                        Console.Write("{0,15:0.0000000}", matrix[i, j]);
                    else
                        Console.Write("{0,25}", matrix[i, j]);
                Console.WriteLine();
            }
            Console.WriteLine();
        }

        /// <summary>
        /// Вывод вектора
        /// </summary>
        /// <param name="matrix"></param>
        /// <param name="n"></param>
        /// <param name="rounded"></param>
        public static void WriteVector(double[] vector, int n, bool rounded = true)
        {
            for (int i = 0; i < n; i++)
            {
                if (!rounded)
                    Console.Write("{0,15:0.0000000}", vector[i]);
                else
                    Console.Write("{0,25}", vector[i]);
            }
            Console.WriteLine();
            Console.WriteLine();
        }

        /// <summary>
        ///  Чтение матрицы из файла
        /// </summary>
        /// <param name="matrix"></param>
        /// <param name="n"></param>
        public static void ReadMatrix(double[,] matrix, int n, string FileName)
        {
            string s;
            int j = 0;
            StreamReader Fr = null;
            try
            {
                Fr = new StreamReader(FileName);
            }
            catch (Exception e)
            {
                Console.WriteLine("Exception:" + e.Message);
            }
            //записываем матрицу в двумерный массив
            while ((s = Fr.ReadLine()) != null)
            {
                for (int i = 0; i < n; i++)
                    matrix[j, i] = double.Parse(s.Split(' ')[i]);
                j++;
            }
            Fr.Close();
        }

        /// <summary>
        /// Умножение матрицы на вектор
        /// </summary>
        /// <param name="A"></param>
        /// <param name="x"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        public static double[] MultMatrixVector(double[,] A, double[] x, int n)
        {
            double[] res = new double[n];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    res[i] += A[i, j] * x[j];
            return res;
        }

        /// <summary>
        /// Умножение вектора на вектор
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        public static double MultVectors(double[] a, double[] b, int n)
        {
            double res = 0;

            for (int i = 0; i < n; i++)
                res += a[i] * b[i];

            return res;
        }

        /// <summary>
        /// Умножение вектора на число
        /// </summary>
        /// <param name="x"></param>
        /// <param name="vector"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        public static double[] MultDigitVector(double x, double[] vector, int n)
        {
            double[] res = new double[n];

            for (int i = 0; i < n; i++)
                res[i] = x * vector[i];

            return res;
        }

        /// <summary>
        /// Проверка диагональных элементов на наличие нуля
        /// </summary>
        /// <param name="A"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        public static double MainDiagonalElements(double[,] A, int n)
        {
            double res = 1;
            for (int i = 0; i < n; i++)
                res *= A[i, i];
            return res;
        }

        /// <summary>
        /// Сложение векторов
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        public static double[] SumVector(double[] A, double[] B, int n)
        {
            for (int i = 0; i < n; i++)
                A[i] += B[i];

            return A;
        }

        /// <summary>
        /// Расчетная формула. Метод простых итераций
        /// </summary>
        /// <param name="args"></param>
        /// 
        public static double[] SimpleIteration(double[,] A, double[] b, double[] PreviousValues, int n)
        {
            double y = 0.9;
            double[] CurrentValues = new double[n];
            if (MainDiagonalElements(A, n) != 0)
            {
                double sum = 0;

                for (int i = 0; i < n; i++)
                {
                    sum = 0;
                    for (int j = 0; j < n; j++)
                    {
                        sum += A[i, j] * PreviousValues[j];
                    }
                    CurrentValues[i] = PreviousValues[i] + (y * 2) / NormMatrix(A, n) * (b[i] - sum);
                }
            }
            return CurrentValues;
        }

        /// <summary>
        /// Расчетная формула. Метод смежных градиентов
        /// </summary>
        /// <param name="args"></param>
        /// 
        public static double[] Gradient(double[,] A, double[] PreviousValues, double[] PreviousPreviousValues,
            double Alpha, double[] b, int n)
        {
            double[] PreviousR = GradientR(A, PreviousPreviousValues, b, n);
            return
                SumVector(
                    SumVector(
                        MultDigitVector(Alpha, PreviousValues, n),
                        MultDigitVector((1 - Alpha), PreviousPreviousValues, n), n),
                    MultDigitVector(
                        -1,
                        MultDigitVector(
                            GradientT(
                                A,
                                GradientR(A, PreviousValues, b, n),
                                n) * Alpha,
                            GradientR(A, PreviousValues, b, n),
                            n),
                        n),
                    n);
        }

        /// <summary>
        /// Расчетная формула. Метод наискорейшего спуска
        /// </summary>
        /// <param name="args"></param>
        /// 
        public static double[] SteepestDescent(double[,] A, double[] PreviousValues, double[] PreviousPreviousValues, double[] b, int n)
        {
            double[] PreviousR = GradientR(A, PreviousPreviousValues, b, n);
            return SumVector(PreviousValues, MultDigitVector(SteepestDescentAlpha(A, PreviousR, n), PreviousR, n), n);
        }

        /// <summary>
        /// Вектор невязки
        /// </summary>
        /// <param name="A"></param>
        /// <param name="PreviousValues"></param>
        /// <param name="b"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        public static double[] GradientR(double[,] A, double[] PreviousValues, double[] b, int n)
        {
            double[] r = MultMatrixVector(A, PreviousValues, n);

            for (int i = 0; i < n; i++)
                r[i] -= b[i];

            return r;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="r"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        public static double GradientT(double[,] A, double[] r, int n)
        {
            return MultVectors(r, r, n) / MultVectors(MultMatrixVector(A, r, n), r, n);
        }

        /// <summary>
        /// Альфа 
        /// </summary>
        /// <param name="PreviousR"></param>
        /// <param name="PreviousPreviousR"></param>
        /// <param name="T"></param>
        /// <param name="PreviousT"></param>
        /// <param name="PreviousAlpha"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        public static double GradientAlpha(double[] PreviousR, double[] PreviousPreviousR, double T,
            double PreviousT, double PreviousAlpha, int n)
        {
            return 1 / (1 - (T * MultVectors(PreviousR, PreviousR, n)) /
                (PreviousT * PreviousAlpha * MultVectors(PreviousPreviousR, PreviousPreviousR, n)));
        }

        /// <summary>
        /// Норма вектора
        /// </summary>
        /// <param name="A"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        public static double NormVector(double[] A, int n)
        {
            double res = 0;

            for (int i = 0; i < n; i++)
                res += A[i] * A[i];

            return Math.Sqrt(res);
        }

        /// <summary>
        /// Оценка нормы матрицы перехода
        /// </summary>
        /// <param name="v"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        public static double AssessmentOfNorms(double[] CurrentValues, double[] PreviousValues, double[] PrePreviousValues, int n)
        {
            double[] A = new double[n];
            double[] B = new double[n];

            for (int i = 0; i < n; i++)
            {
                A[i] = Math.Abs(CurrentValues[i] - PreviousValues[i]);
                B[i] = Math.Abs(PreviousValues[i] - PrePreviousValues[i]);
            }

            return NormVector(A, n) / NormVector(B, n);
        }

        /// <summary>
        /// Копирование вектора
        /// </summary>
        /// <param name="M"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        public static double[] CopyVector(double[] v, int n)
        {
            double[] res = new double[n];
            for (int i = 0; i < n; i++)
                res[i] = v[i];
            return res;
        }

        /// <summary>
        /// Умножение матриц
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        public static double[,] MultMatrix(double[,] A, double[,] B, int n)
        {
            double[,] AB = new double[n, n];
            for (int t = 0; t < n; t++)
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++)
                        AB[t, i] += A[t, j] * B[j, i];
            return AB;
        }

        /// <summary>
        /// Транспонированная матрица
        /// </summary>
        /// <param name="A"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        public static double[,] RevertMatrix(double[,] A, int n)
        {
            double[,] res = new double[n, n];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    res[i, j] = A[j, i];
            return res;
        }

        /// <summary>
        /// Норма матрицы
        /// </summary>
        /// <param name="A"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        public static double NormMatrix(double[,] A, int n)
        {
            double maxlambda = double.MinValue;
            double[] lambda;
            double[,] z;
            alglib.smatrixevd(MultMatrix(RevertMatrix(A, n), A, n), n, 0, false, out lambda, out z);
            for (int i = 0; i < n; i++)
                maxlambda = Math.Max(maxlambda, Math.Abs(lambda[i]));
            return Math.Sqrt(maxlambda);
        }        
        
        /// <summary>
        /// Альфа 
        /// </summary>
        /// <param name="PreviousR"></param>
        /// <param name="PreviousPreviousR"></param>
        /// <param name="T"></param>
        /// <param name="PreviousT"></param>
        /// <param name="PreviousAlpha"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        public static double SteepestDescentAlpha(double[,] A, double[] PreviousR, int n)
        {
            return MultVectors(PreviousR, PreviousR, n) / MultVectors(MultMatrixVector(A, PreviousR, n), PreviousR, n);
        }

        /// <summary>
        /// Метод простых итераций
        /// </summary>
        /// <param name="A"></param>
        /// <param name="b"></param>
        /// <param name="x"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        public static double[] MetodOfSimpleIteration(double[,] A, double[] b, int n, double eps)
        {
            double y = 0.9;
            double[] PrePreviousValues = new double[n];
            double[] PreviousValues = new double[n];
            double[] CurrentValues = new double[n];
            int k = 1;
            while (true)
            {
                // значения неизвестных на текущей итерации
                CurrentValues = SimpleIteration(A, b, PreviousValues, n);

                // текущая погрешность относительно предыдущей итерации
                double delta = AssessmentOfNorms(CurrentValues, PreviousValues, PrePreviousValues, n);

                WriteInformation(k, A, PreviousValues, b, n, delta);
                Console.WriteLine("Тета = " + (y * 2) / NormMatrix(A, n));
                Console.WriteLine();

                if (delta < eps)
                    break;
                //Переходим к следующей итерации. Предыдущие значения неизвестных становятся значениями на предпредыдущей итерации
                PrePreviousValues = CopyVector(PreviousValues, n);
                //Текущие значения неизвестных становятся значениями на предыдущей итерации
                PreviousValues = CopyVector(CurrentValues, n);

                k++;
            }
            Console.WriteLine("Количество итераций = " + k);
            Console.WriteLine();
            return PreviousValues;
        }

        /// <summary>
        /// Метод сопряженных градиентов
        /// </summary>
        /// <param name="A"></param>
        /// <param name="b"></param>
        /// <param name="n"></param>
        /// <param name="eps"></param>
        /// <returns></returns>
        public static double[] MetodGradient(double[,] A, double[] b, int n, double eps)
        {
            double[] PrePreviousValues = new double[n];
            double[] PreviousValues = new double[n];
            double[] CurrentValues = new double[n];
            int k = 1;
            double Alpha = 1;
            while (true)
            {
                // значения неизвестных на текущей итерации
                CurrentValues = Gradient(A, PreviousValues, PrePreviousValues, Alpha, b, n);

                // текущая погрешность относительно предыдущей итерации
                double delta = AssessmentOfNorms(CurrentValues, PreviousValues, PrePreviousValues, n);

                WriteInformation(k, A, PreviousValues, b, n, delta);
                Console.WriteLine("Альфа = " + Alpha);
                Console.WriteLine();

                if (delta < eps)
                    break;
                //Переходим к следующей итерации. Предыдущие значения неизвестных становятся значениями на предпредыдущей итерации
                PrePreviousValues = CopyVector(PreviousValues, n);
                //Текущие значения неизвестных становятся значениями на предыдущей итерации
                PreviousValues = CopyVector(CurrentValues, n);

                double[] R = GradientR(A, PreviousValues, b, n);
                double[] PreviousR = GradientR(A, PrePreviousValues, b, n);

                Alpha = GradientAlpha(
                    R,
                    PreviousR,
                    GradientT(A, R, n),
                    GradientT(A, PreviousR, n),
                    Alpha,
                    n);

                k++;
            }
            Console.WriteLine("Количество итераций = " + k);
            Console.WriteLine();
            return PreviousValues;
        }

        /// <summary>
        /// Метод наискорейшего спуска
        /// </summary>
        /// <param name="args"></param>
        public static double[] MetodSteepestDescent(double[,] A, double[] b, int n, double eps)
        {
            double[] PrePreviousValues = new double[n];
            double[] PreviousValues = new double[n];
            double[] CurrentValues = new double[n];
            int k = 1;
            while (true)
            {
                // значения неизвестных на текущей итерации
                CurrentValues = SteepestDescent(A, PreviousValues, PrePreviousValues, b, n);

                // текущая погрешность относительно предыдущей итерации
                double delta = AssessmentOfNorms(CurrentValues, PreviousValues, PrePreviousValues, n);

                WriteInformation(k, A, PreviousValues, b, n, delta);

                if (delta < eps)
                    break;
                //Переходим к следующей итерации. Предыдущие значения неизвестных становятся значениями на предпредыдущей итерации
                PrePreviousValues = CopyVector(PreviousValues, n);
                //Текущие значения неизвестных становятся значениями на предыдущей итерации
                PreviousValues = CopyVector(CurrentValues, n);

                double[] R = GradientR(A, PreviousValues, b, n);
                double[] PreviousR = GradientR(A, PrePreviousValues, b, n);  
            }
            Console.WriteLine("Количество итераций = " + k);
            Console.WriteLine();
            return PreviousValues;
        }

        static void Main(string[] args)
        {
            double eps = 0.0001;
            int n = 4;//размерность матрицы
            double[] x = { 1, 2, 3, 4 }; // вектор x
            double[] b = new double[n]; // вектор b
            double[,] A = new double[n, n];//матрица А
            double[] SI = new double[n]; // вектор x
            double[] G = new double[n]; // вектор x
            double[] SD = new double[n]; // вектор x
            double[] PVR = new double[n]; // вектор x

            ReadMatrix(A, n, "Input10.txt");// считываем матрицу A из файл

            Console.WriteLine("eps:" + eps);// выведем e

            Console.WriteLine("x:");// выведем вектор x
            WriteVector(x, n);

            Console.WriteLine("A: ");// выведем матрицу А
            WriteMatrix(A, n);

            b = MultMatrixVector(A, x, n); //найдем решение уравнения
            Console.WriteLine("b: ");// выведем вектор b
            WriteVector(b, n);

            Console.WriteLine("Метод простых итераций: ");
            SI = MetodOfSimpleIteration(A, b, n, eps);
            Console.WriteLine("x: ");// выведем вектор X
            WriteVector(SI, n);

            Console.WriteLine("Метод сопряженных градиентов: ");
            G = MetodGradient(A, b, n, eps);
            Console.WriteLine("x: ");// выведем вектор X
            WriteVector(G, n);

            //Console.WriteLine("Метод наискорейшего спуска: ");
            //SD = MetodSteepestDescent(A, b, n, eps);
            //Console.WriteLine("x: ");// выведем вектор X
            //WriteVector(SD, n);

            Console.ReadKey();
        }
    }
}