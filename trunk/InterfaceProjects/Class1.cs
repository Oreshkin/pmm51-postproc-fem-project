using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace fem_interface
{
    struct point
    {
        int x;
        int y;
    }
    /// <summary>
    /// N- размерность матрицы
    /// </summary>
    struct Imatrix
    {
        int N;
        double[,] Matrix;
    };
    public class Class1
    {
    }
    public interface IBasisMKE
    {
        /// <summary>
        /// 
        /// </summary>
        /// <param name="gr"></param>
        /// <param name="number" - номер элемента сетки></param>
        void drawElem(IGrafic gr, int number);
        Imatrix LocalMatrix(int matrixType, int BasisType);
        double GetValuePoint(point A, int BasisType);
    }
    class triangle : IBasisMKE
    {
        /// <summary>
        /// значение базисной функции в точке
        /// </summary>
        /// <param name="number" - номер базисной функции></param>
        /// <param name="A" - точка, в которой нужно значение></param>
        /// <param name="BasisType" - тип базиса (линейный,  квадрат., кубич.)></param>
        /// <returns></returns>
        double basis(int number, point A, int BasisType)
        {
            double q;
            switch (BasisType)
            {
                case (1):
                    {
                        switch (number)
                        {//линейный
                            case (1): { q= break; }
                            case (2): {  break; }
                            case (3): {  break; }
                            default: { break; }
                        }
                        break;
                    }
                case (2):
                    {   // квадратичный
                        switch (number)
                        {
                            case (1): { A = Mass_sqr(); break; }
                            case (2): { A = Gest_sqr(); break; }
                            case (3): { A = Exotic_sqr(); break; }
                            default: { break; }
                        }
                        break;
                    }
                case (3):
                    {   // кубический 
                        switch (number)
                        {
                            case (1): { A = Mass_cub(); break; }
                            case (2): { A = Gest_cub(); break; }
                            case (3): { A = Exotic_cub(); break; }
                            default: { break; }
                        }
                        break;
                    }
                default: { break; }
            }
            return q;
        }
        public void drawElem(IGrafic gr, int number)
        { 
        
        }
        public Imatrix LocalMatrix(int matrixType, int BasisType)
        {
            Imatrix A;
            switch (BasisType)
            {
                case (1): 
                {
                    switch (matrixType)
                    {
                        case (1): { A = Mass_lin(); break; }
                        case (2): { A = Gest_lin(); break; }
                        case (3): { A = Exotic_lin(); break; }
                        default: { break;}
                    }
                    break;
                }
                case (2): 
                    {
                        switch (matrixType)
                        {
                            case (1): { A = Mass_sqr(); break; }
                            case (2): { A = Gest_sqr(); break; }
                            case (3): { A = Exotic_sqr(); break; }
                            default: { break; }
                        }
                        break;
                    }
                case (3): 
                    {
                        switch (matrixType)
                        {
                            case (1): { A = Mass_cub(); break; }
                            case (2): { A = Gest_cub(); break; }
                            case (3): { A = Exotic_cub(); break; }
                            default: { break; }
                        }
                        break;
                    }
                default: { break; }
            }
            return A;
        }
        public double GetValuePoint(point A, int BasisType)
        { 
        
        }
    
    }
    class rectangle : IBasisMKE
    {
        /// <summary>
        /// значение базисной функции в точке
        /// </summary>
        /// <param name="number" - номер базисной функции></param>
        /// <param name="A" - точка, в которой нужно значение></param>
        /// <param name="BasisType" - тип базиса (линейный,  квадрат., кубич.)></param>
        /// <returns></returns>
        double basis(int number, point A, int BasisType)
        {
            double q;
            switch (BasisType)
            {
                case (1):
                    {
                        switch (number)
                        {//линейный
                            case (1): { q= break; }
                            case (2): {  break; }
                            case (3): {  break; }
                            default: { break; }
                        }
                        break;
                    }
                case (2):
                    {   // квадратичный
                        switch (number)
                        {
                            case (1): { A = Mass_sqr(); break; }
                            case (2): { A = Gest_sqr(); break; }
                            case (3): { A = Exotic_sqr(); break; }
                            default: { break; }
                        }
                        break;
                    }
                case (3):
                    {   // кубический 
                        switch (number)
                        {
                            case (1): { A = Mass_cub(); break; }
                            case (2): { A = Gest_cub(); break; }
                            case (3): { A = Exotic_cub(); break; }
                            default: { break; }
                        }
                        break;
                    }
                default: { break; }
            }
            return q;
        }
        public void drawElem(IGrafic gr, int number)
        {

        }
        public Imatrix LocalMatrix(int matrixType, int BasisType)
        {
            Imatrix A;
            switch (BasisType)
            {
                case (1):
                    {
                        switch (matrixType)
                        {
                            case (1): { A = Mass_lin(); break; }
                            case (2): { A = Gest_lin(); break; }
                            case (3): { A = Exotic_lin(); break; }
                            default: { break; }
                        }
                        break;
                    }
                case (2):
                    {
                        switch (matrixType)
                        {
                            case (1): { A = Mass_sqr(); break; }
                            case (2): { A = Gest_sqr(); break; }
                            case (3): { A = Exotic_sqr(); break; }
                            default: { break; }
                        }
                        break;
                    }
                case (3):
                    {
                        switch (matrixType)
                        {
                            case (1): { A = Mass_cub(); break; }
                            case (2): { A = Gest_cub(); break; }
                            case (3): { A = Exotic_cub(); break; }
                            default: { break; }
                        }
                        break;
                    }
                default: { break; }
            }
            return A;
        }
        public double GetValuePoint(point A, int BasisType)
        {

        }

    }
   
}
