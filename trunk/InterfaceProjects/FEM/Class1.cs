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
        void GetNodes(point[] p, int num); // получение КЭ по глобальному номеру 
        MaterialIdentifire GetMaterial(int num); // получение материала по глобальному номеру КЭ 
        double GetValues(point A); //получение решения  в точке A 
        double[,] LocalMatrix(MatrixType type, int num); // получение локальной матрицы заданного типа type для элемента  с глобальным номером num 
    }
    public class triangleLin : IBasisMKE  //линейные треугольники 
    {
        public Mesh Set;     // сетка считывается и создается в конструкторе объекта этого класса
        public triangleLin(string FileNameMesh, string FileNameSolution)
        {
            Set = new Mesh(FileNameMesh,FileNameSolution);//чтение сетки
        }
        public void GetNodes(point[] p, int num)   // получили по номеру элемента координаты вершин
        {
            int i;
            for (i = 0; i < 3; i++)
            {
                p[i].x = Set.elements[num].vertex[i].x;
                p[i].y = Set.elements[num].vertex[i].y;
            }
        }
        public MaterialIdentifire GetMaterial(int num)   // получение номера материала
        { 
            return Set.elements[num].material; 
        }
        public double GetValues(point A)  // получение решения в точке A
        {
            int i;
            int num=-1;//номер элемента
            double S1, S2, S3, S;
            double a1, a2, a3, b1, b2, b3, c1, c2, c3;
            double x=A.x;
            double y=A.y;
            double x1,x2,x3,y1,y2,y3;
            double[] w = new double[3];  // базисные функции
            //поиск нужного конечного элемента
            for (i = 0; i < Set.elements.Length; i++)
            {
                x1=Set.elements[i].vertex[0].x;
                x2=Set.elements[i].vertex[1].x;
                x3=Set.elements[i].vertex[2].x;
                y1=Set.elements[i].vertex[0].y;
                y2=Set.elements[i].vertex[1].y;
                y3=Set.elements[i].vertex[2].y;
                S=Math.Abs((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1));
                S1 = Math.Abs((x2 - x) * (y3 - y) - (x3 - x) * (y2 - y));
                S2 = Math.Abs((x - x1) * (y3 - y1) - (x3 - x1) * (y - y1));
                S3 = Math.Abs((x2 - x1) * (y - y1) - (x - x1) * (y2 - y1));
                if (Math.Abs(S1 + S2 + S3 - S) < 1e-10)
                {
                    num = i;  // нашли
                    break;
                }
            }
            if (num == -1) return 0;//ошибка, не нашли элемент, которому принадлежит точка А
            else
            {
                x1 = Set.elements[num].vertex[0].x;
                x2 = Set.elements[num].vertex[1].x;
                x3 = Set.elements[num].vertex[2].x;
                y1 = Set.elements[num].vertex[0].y;
                y2 = Set.elements[num].vertex[1].y;
                y3 = Set.elements[num].vertex[2].y;
                S = Math.Abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));
                //барицентрические координаты
                a1 = (x2 * y3 - x3 * y2)/S;
                a2 = (x3 * y1 - x1 * y3)/S;
                a3 = (x1 * y2 - x2 * y1)/S;

                b1 = (y2 - y3)/S;
                b2 = (y3 - y1)/S;
                b3 = (y1 - y2)/S;

                c1 = (x3 - x2)/S;
                c2 = (x1 - x3)/S;
                c3 = (x2 - x1)/S;
                //барицентрические координаты
                //базисные функции в точке
                w[0] = a1 + b1 * x + c1 * y;
                w[1] = a2 + b2 * x + c2 * y;
                w[2] = a3 + b3 * x + c3 * y;

                return w[0] * Set.solution[Set.elements[num].vertex[0].globalNum] + w[1] * Set.solution[Set.elements[num].vertex[1].globalNum] + w[2] * Set.solution[Set.elements[num].vertex[2].globalNum];
            }

        }
        public double[,] LocalMatrix(MatrixType type, int num)
        {
            double[,]M=new double[4,4];
            int i, j;
            double x1,x2,x3,y1,y2,y3;
            double det;//определитель
            double gamma = Set.MaterialGamma(Set.elements[num].material);
            if (type == MatrixType.mass)
            {
                x1=Set.elements[num].vertex[0].x;
                x2=Set.elements[num].vertex[1].x;
                x3=Set.elements[num].vertex[2].x;
                y1=Set.elements[num].vertex[0].y;
                y2=Set.elements[num].vertex[1].y;
                y3=Set.elements[num].vertex[2].y;
                det=Math.Abs((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1));
                for (i = 0; i < 3; i++)
                    for (j = 0; j < 3; j++)
                        if (i == j)
                            M[i,j] = gamma*det/12.0;
                        else
                            M[i,j] = gamma*det/24.0;
            }
            else
            {
                   
            }
           return M;
        }// получение локальной матрицы заданного типа type для элемента  с глобальным номером num 
    }

   /* public class triangleQuard : IBasisMKE //квадратичные треугольники 
    {
        Mesh Set;
        public void GetNodes(point[] p, int num)
        {
            int i;
            for (i = 0; i < 3; i++)
            {
                p[i].x = Set.elements[num].vertex[i].x;
                p[i].y = Set.elements[num].vertex[i].y;
            }
        }
        public int GetMaterial(int num)
        { return 1; }
        public void GetValues(double[] q, point A, int num)
        { }
        public double[,] LocalMatrix(MatrixType type, int num)
        {
            double[,]M=new double[4,4];
            return M;
        }// получение локальной матрицы заданного типа type для элемента  с глобальным номером num 
    }

    public class triangleCub : IBasisMKE //кубические треугольники 
    {
        Mesh Set;
        public void GetNodes(point[] p, int num)
        {
            int i;
            for (i = 0; i < 3; i++)
            {
                p[i].x = Set.elements[num].vertex[i].x;
                p[i].y = Set.elements[num].vertex[i].y;
            }
        }
        public int GetMaterial(int num)
        { return 1; }
        public void GetValues(double[] q, point A, int num)
        { }
        public double[,] LocalMatrix(MatrixType type, int num)
        {
            double[,]M=new double[4,4];
            return M;
        }// получение локальной матрицы заданного типа type для элемента  с глобальным номером num 
    }*/
    public struct point
    {
        public double x;
        public double y;
        public int globalNum; // использовать только для вершин сетки/// глобальный номер узла в сетке
    }

    public enum ElemType { triangle, rectangle };//тип элемента:{треугольник, четырехугольник} 
    public enum BasisType { lin, quadr, cub }; // тип базиса:{линейный, квадратичный, тубический} 
    public enum MatrixType { mass, gest, exotic };//тип матрицы:{массы, жесткости, особые} 
    public enum MaterialIdentifire { one, two, three };//номер подобласти/материала 


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
        public Element[] elements;  //массив элементов 
        public double[] solution;   // массив весов 
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
                    elements[i].vertex[j].globalNum = qwe;
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
        public double MaterialGamma( MaterialIdentifire type )
        {
            if(type==MaterialIdentifire.one) 
                return 1.0;
            else
                if(type==MaterialIdentifire.two)
                    return 2.0;
                else
                    return 3.0;
        }
    } 
    
} 




