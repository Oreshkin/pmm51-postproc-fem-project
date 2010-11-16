namespace FEM
{
using System;
using System.Collections.Generic;
using System.Collections;
using System.Linq;
using System.Text;
using System.IO;

    public interface IBasisMKE
    {
        void GetNodes(point[] p, int num); // полечение КЭ по глобальному номеру
        int GetMaterial(int num); // получение материала по глобальному номеру КЭ
        void GetValues(double[] q, int num); //получение значения 
    }
    public class triangleLin : IBasisMKE  //линейные треугольники
    {
        public void GetNodes(point[] p, int num)
        { }
        public int GetMaterial(int num)
        { return 1; }
        public void GetValues(double[] q, int num)
        { }
    }

    public class triangleQuard : IBasisMKE //квадратичные треугольники
    {
        public void GetNodes(point[] p, int num)
        { }
        public int GetMaterial(int num)
        { return 1; }
        public void GetValues(double[] q, int num)
        { }
    }

    public class triangleCub : IBasisMKE //кубические треугольники
    {
        public void GetNodes(point[] p, int num)
        { }
        public int GetMaterial(int num)
        { return 1; }
        public void GetValues(double[] q, int num)
        { }
    }
    public struct point
    {
        public double x;
        public double y;
    }

    public enum ElemType { triangle, rectangle };//тип элемента:{треугольник, четырехугольник}
    public enum BasisType { lin, quadr, cub }; // тип базиса:{линейный, квадратичный, тубический}
    public enum MatrixType { mass, gest, exotic };//тип матрицы:{массы, жесткости, особые}
    public enum MaterialIdentifire { one, two,three };//номер подобласти/материала

    public class Element
    {
        public ElemType type;
        public BasisType btype;
        public point[] vertex; // массив вершин КЭ
        public MaterialIdentifire material;

        public Element()
        {

        }

        public void ElementCreate()
        {
            if (type == ElemType.triangle)

                vertex = new point[3];

            else
                vertex = new point[4];
        }

    }


    public class Mesh : IEnumerable, IEnumerator
    {
        Element[] elements;  //массив элементов
        double[] solution;   // массив весов
        public Mesh(string FileNameMesh, string FileNameSolution)  //чтение сетки и решения из файлов
        {
            FileStream f = new FileStream(FileNameMesh, FileMode.Open, FileAccess.Read);
            StreamReader Reader = new StreamReader(f);
            string[] t;
            int i, N, N_elem, qwe, flag_of_vertex, j;
            string[] temp;

            for (i = 0; i < 4; i++)
                t = Reader.ReadLine().Split(' ');
            t = Reader.ReadLine().Split(' ');
            N = Convert.ToInt32(t[0]);
            point[] nodes = new point[N];
            for (i = 0; i < N; i++)
            {
                temp = Reader.ReadLine().Split(' ');
                nodes[i].x = Convert.ToDouble(temp[1]);
                nodes[i].y = Convert.ToDouble(temp[2]);
            }
            t = Reader.ReadLine().Split(' ');
            t = Reader.ReadLine().Split(' ');
            N_elem = Convert.ToInt32(t[0]);
            elements = new Element[N_elem];
            for (i = 0; i < N_elem; i++)
            {
                temp = Reader.ReadLine().Split(' ');
                if (Convert.ToInt32(temp[1]) == 2)
                {
                    elements[i].type = ElemType.triangle;
                    flag_of_vertex = 3;
                }
                else
                {
                    elements[i].type = ElemType.rectangle;
                    flag_of_vertex = 4;
                }
                if (Convert.ToInt32(temp[3]) == 1)
                    elements[i].material = MaterialIdentifire.one;
                else
                    if (Convert.ToInt32(temp[3]) == 2)
                        elements[i].material = MaterialIdentifire.two;
                    else
                        elements[i].material = MaterialIdentifire.three;
                elements[i].ElementCreate();
                for (j = 0; j < flag_of_vertex; j++)
                {
                    qwe = Convert.ToInt32(temp[6 + j]);
                    elements[i].vertex[j].x = nodes[qwe].x;
                    elements[i].vertex[j].y = nodes[qwe].y;
                }

            }

            Reader.Close();

            solution = new double[N];
            f = new FileStream(FileNameSolution, FileMode.Open, FileAccess.Read);
            Reader = new StreamReader(f);
            for (i = 0; i < N; i++)
            {
                temp = Reader.ReadLine().Split(' ');
                solution[i] = Convert.ToDouble(temp[0]);

            }
            Reader.Close();
        }

        public double solution_max()
        {
            double max = solution[0];
            for (int i = 1; i < solution.Length; i++)
                if (max < solution[i]) max = solution[i];
            return max;
        }

        public double solution_min()
        {
            double min = solution[0];
            for (int i = 1; i < solution.Length; i++)
                if (min > solution[i]) min = solution[i];
            return min;
        }
        public IEnumerator GetEnumerator()
        {
            return (IEnumerator)this;
        }
        int pos = -1;

        public bool MoveNext()
        {
            if (pos < elements.Length - 1)
            {
                pos++;
                return true;
            }
            else
                return false;
        }
        public object Current
        {
            get { return elements[pos]; }
        }
        public void Reset()
        {
            pos = -1;
        }
    }
}
