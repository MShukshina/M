using System;
using System.IO;
using System.Threading.Tasks;


namespace Lab3_Чисмет
{
    class Program
    {
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
        }

        /// <summary>
        ///  Чтение матрицы из файла
        /// </summary>
        /// <param name="matrix"></param>
        /// <param name="n"></param>
        public static void ReadMatrix(double[,] matrix, int n)
        {
            string s;
            int j = 0;
            StreamReader Fr = null;
            try
            {
                Fr = new StreamReader("Input10.txt");
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
        /// Проверка диагональных элементов на наличие нуля
        /// </summary>
        /// <param name="A"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        public static double MainDiagonalElements(double[,] A, int n)
        {
            double res = 0;
            for (int i = 0; i < n; i++)
                    res += A[i, i];
            return res;
        }

        /// <summary>
        /// Приведим систему к виду, удобному для итерации
        /// </summary>
        /// <param name="args"></param>
        /// 
        public static double[] SimpleIteration(double[,] A, double[] b,double[] PreviousValues, int n)
        {
            double y = 0.9;
            double[] CurrentValues = new double[n];
            if (MainDiagonalElements(A,n) != 0)
            {
                double sum = 0;
                
                for (int i = 0; i < n; i++)
                {
                    sum = 0;
                    for (int j = 0; j < n; j++)
                    {
                            sum += A[i, j] * PreviousValues[j];
                    }
                    CurrentValues[i] = PreviousValues[i]  +  (y*2)/ NormMatrix(A,n) * (b[i] - sum);
                }
            }
            return CurrentValues;
        }

        /// <summary>
        /// Оценка нормы матрицы перехода
        /// </summary>
        /// <param name="v"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        public static double AssessmentOfNorms(double[] CurrentValues, double[] PreviousValues,double[] PrePreviousValues, int n)
        {
            double q = 0;
            for (int i = 0; i < n; i++)
                q =(Math.Abs(CurrentValues[i] - PreviousValues[i])) / Math.Abs(PreviousValues[i] - PrePreviousValues[i]);
            return q;
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
        /// Метод простых итераций
        /// </summary>
        /// <param name="A"></param>
        /// <param name="b"></param>
        /// <param name="x"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        public static double[] MetodOfSimpleIteration(double[,] A, double[] b, int n, double eps)
        {
            double[] PrePreviousValues = new double[n];
            double[] PreviousValues = new double[n];
            double[] CurrentValues = new double[n];
            while (true)
            {
                // значения неизвестных на текущей итерации
                CurrentValues = SimpleIteration(A,b,PreviousValues,n);

                // текущая погрешность относительно предыдущей итерации
                double error = AssessmentOfNorms(CurrentValues, PreviousValues, PrePreviousValues, n);

                if (error < eps)
                        break;
                //Переходим к следующей итерации. Предыдущие значения неизвестных становятся значениями на предпредыдущей итерации
                PrePreviousValues = CopyVector(PreviousValues,n);
                //Текущие значения неизвестных становятся значениями на предыдущей итерации
                PreviousValues = CopyVector(CurrentValues, n);
            }
            return PreviousValues;
        }

            static void Main(string[] args)
        {
            double eps = 0.0001;
            int n = 4;//размерность матрицы
            double[] x = { 1, 2, 3, 4 }; // вектор x
            double[] b = new double[n]; // вектор b
            double[,] A = new double[n,n];//матрица А
            double[] X = new double[n]; // вектор x. Метод Якоби


            ReadMatrix(A, n);// считываем матрицу A из файл

            Console.WriteLine("eps:" + eps);// выведем e

            Console.WriteLine("x:");// выведем вектор x
            WriteVector(x, n);

            Console.WriteLine("A: ");// выведем матрицу А
            WriteMatrix(A, n);

            b = MultMatrixVector(A, x, n); //найдем решение уравнения
            Console.WriteLine("b: ");// выведем вектор b
            WriteVector(b, n);

            X = MetodOfSimpleIteration(A, b, n, eps);
            Console.WriteLine("X: ");// выведем вектор X
            WriteVector(X, n);

            Console.ReadKey();
        }
    }
}
