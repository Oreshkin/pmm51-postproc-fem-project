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
        MaterialIdentifire GetMaterial(int num); // получение материала по глобальному номеру КЭ 
        double GetValues(point A, TypeOfSolution type, Direction d); //получение решения  в точке A 
        double[,] LocalMatrix(MatrixType type, int num); // получение локальной матрицы заданного типа type для элемента  с глобальным номером num 
        double[] DiffDirection(Direction d); //смена направления производной, d-переменная, по которой дифференцируем d={x,y}
        double solution_max(TypeOfSolution type, Direction d); //поиск максимума
        double solution_min(TypeOfSolution type, Direction d); //поиск минимума
    }
    public class triangleLin : IBasisMKE  //линейные треугольники 
    {
        public int flagX, flagY;
        public Mesh Set;     // сетка считывается и создается в конструкторе объекта этого класса
        public SLAE_MSG S_X; // решатель для dx
        public SLAE_MSG S_Y; // решатель для dy
        public triangleLin(string FileNameMesh, string FileNameSolution, string FileNameIG, string FileNameJG)
        {
            Set = new Mesh(FileNameMesh, FileNameSolution);//чтение сетки
            S_X = new SLAE_MSG(FileNameIG, FileNameJG, Set.solution.Length);
            S_Y = new SLAE_MSG(FileNameIG, FileNameJG, Set.solution.Length);
        }
        public void GetNodes(point[] p, int num)   // получили по номеру элемента координаты вершин
        {
            int i;
            for (i = 0; i < p.Length; i++) // было до 3-х
            {
                p[i].x = Set.elements[num].vertex[i].x;
                p[i].y = Set.elements[num].vertex[i].y;
                p[i].globalNum = Set.elements[num].vertex[i].globalNum;
            }
        }
        public MaterialIdentifire GetMaterial(int num)   // получение номера материала
        {
            return Set.elements[num].material;
        }
        public double GetValues(point A, TypeOfSolution type, Direction d)  // получение решения в точке A
        {
            int i;
            int num = -1;//номер элемента
            double S1, S2, S3, S;
            double a1, a2, a3, b1, b2, b3, c1, c2, c3;
            double x = A.x;
            double y = A.y;
            double x1, x2, x3, y1, y2, y3;
            double[] w = new double[3];  // базисные функции
            //поиск нужного конечного элемента
            for (i = 0; i < Set.elements.Length; i++)
            {
                x1 = Set.elements[i].vertex[0].x;
                x2 = Set.elements[i].vertex[1].x;
                x3 = Set.elements[i].vertex[2].x;
                y1 = Set.elements[i].vertex[0].y;
                y2 = Set.elements[i].vertex[1].y;
                y3 = Set.elements[i].vertex[2].y;
                S = Math.Abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));
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
                a1 = (x2 * y3 - x3 * y2) / S;
                a2 = (x3 * y1 - x1 * y3) / S;
                a3 = (x1 * y2 - x2 * y1) / S;

                b1 = (y2 - y3) / S;
                b2 = (y3 - y1) / S;
                b3 = (y1 - y2) / S;

                c1 = (x3 - x2) / S;
                c2 = (x1 - x3) / S;
                c3 = (x2 - x1) / S;
                //барицентрические координаты
                //базисные функции в точке
                w[0] = a1 + b1 * x + c1 * y;
                w[1] = a2 + b2 * x + c2 * y;
                w[2] = a3 + b3 * x + c3 * y;
                if (type == TypeOfSolution.Solution)
                    return w[0] * Set.solution[Set.elements[num].vertex[0].globalNum] + w[1] * Set.solution[Set.elements[num].vertex[1].globalNum] + w[2] * Set.solution[Set.elements[num].vertex[2].globalNum];
                else
                {
                    if(d==Direction.x)
                        return w[0] * S_X.Q[Set.elements[num].vertex[0].globalNum] + w[1] * S_X.Q[Set.elements[num].vertex[1].globalNum] + w[2] * S_X.Q[Set.elements[num].vertex[2].globalNum];
                    else return w[0] * S_Y.Q[Set.elements[num].vertex[0].globalNum] + w[1] * S_Y.Q[Set.elements[num].vertex[1].globalNum] + w[2] * S_Y.Q[Set.elements[num].vertex[2].globalNum];
                }
            }

        }
        public double[,] LocalMatrix(MatrixType type, int num)
        {
            double[,] M = new double[3, 3];
            int i, j;
            double[] b = new double[3]; 
            double[] c = new double[3];
            double x1, x2, x3, y1, y2, y3;
            double det;//определитель
            double gamma = Set.MaterialGamma(Set.elements[num].material);

            x1 = Set.elements[num].vertex[0].x;
            x2 = Set.elements[num].vertex[1].x;
            x3 = Set.elements[num].vertex[2].x;
            y1 = Set.elements[num].vertex[0].y;
            y2 = Set.elements[num].vertex[1].y;
            y3 = Set.elements[num].vertex[2].y;
            det = Math.Abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));

            b[0] = (y2 - y3) / det;
            b[1] = (y3 - y1) / det;
            b[2] = (y1 - y2) / det;

            c[0] = (x3 - x2) / det;
            c[1] = (x1 - x3) / det;
            c[2] = (x2 - x1) / det;

            if (type == MatrixType.mass)
            {
                for (i = 0; i < 3; i++)
                    for (j = 0; j < 3; j++)
                        if (i == j)
                            M[i, j] = gamma * det / 12.0;
                        else
                            M[i, j] = gamma * det / 24.0;
            }
            else
            {
                if (type == MatrixType.exotic1)
                {
                    // посчитать матрицы для производной по х
                    for (i = 0; i < 3; i++)
                        for (j = 0; j < 3; j++)
                                M[i, j] = b[j] * det / 6.0;
                }
                else
                {
                    // посчитать матрицы для производной по У
                    for (i = 0; i < 3; i++)
                        for (j = 0; j < 3; j++)
                                M[i, j] = c[j] * det / 6.0;
                }
            }
            return M;
        }// получение локальной матрицы заданного типа type для элемента  с глобальным номером num 

        public double solution_max(TypeOfSolution type, Direction d)
        {
            double max;
            if (type == TypeOfSolution.Solution)
            {
                max = Set.solution[0];
                for (int i = 1; i < Set.solution.Length; i++)
                    if (max < Set.solution[i]) max = Set.solution[i];
                return max;
            }
            else
            {
                if (d == Direction.x)
               {
                    max = S_X.Q[0];
                    for (int i = 1; i < S_X.Q.Length; i++)
                        if (max < S_X.Q[i]) max = S_X.Q[i];
                    return max;
                } 
                else
                {
                    max = S_Y.Q[0];
                    for (int i = 1; i < S_Y.Q.Length; i++)
                        if (max < S_Y.Q[i]) max = S_Y.Q[i];
                    return max;
                } 
            }
        }

        public double solution_min(TypeOfSolution type, Direction d)
        {
            double min;
            if (type == TypeOfSolution.Solution)
            {
                min = Set.solution[0];
                for (int i = 1; i < Set.solution.Length; i++)
                    if (min > Set.solution[i]) min = Set.solution[i];
                return min;
            }
            else
            {
                if (d == Direction.x)
                {
                    min = S_X.Q[0];
                    for (int i = 1; i < S_X.Q.Length; i++)
                        if (min > S_X.Q[i]) min = S_X.Q[i];
                    return min;
                }
                else
                {
                    min = S_Y.Q[0];
                    for (int i = 1; i < S_Y.Q.Length; i++)
                        if (min > S_Y.Q[i]) min = S_Y.Q[i];
                    return min;
                }
            }
        }
        public double[] DiffDirection(Direction d)
        {
            int nkel, i, j;
            double t;
            point[] Velem=new point[3];
            double[,] M = new double[3, 3];
            double[,] Ex = new double[3, 3];
            if (d == Direction.x)
            {
                if (flagX == 0)
                {
                    for(i=0;i<Set.solution.Length;i++) 
                        S_X.F[i]=0.0;

                    for (nkel = 0; nkel < Set.elements.Length; nkel++)  // сборка глобального вектора и матрицы
                    {
                        GetNodes(Velem, nkel);
                        M = LocalMatrix(MatrixType.mass, nkel);
                        Ex = LocalMatrix(MatrixType.exotic1, nkel);
                        SR_GM_INGLOB(Velem, M, Direction.x);
                        for (i = 0; i < 3; i++)
                        {
                            for (j = 0; j < 3; j++)
                            {
                                t = GetValues(Velem[j], TypeOfSolution.Solution, Direction.x);
                                S_X.F[Velem[i].globalNum] += Ex[i, j] * t;
                            }
                        }
                            
                    }
                    S_X.MSG();
                    flagX = 1;
                }
                return S_X.Q;
            }
            else
            {
                if (flagY == 0)
                {
                    for (i = 0; i < Set.solution.Length; i++)
                        S_Y.F[i] = 0.0;

                    for (nkel = 0; nkel < Set.elements.Length; nkel++) // сборка глобального вектора и матрицы
                    {
                        GetNodes(Velem, nkel);
                        M = LocalMatrix(MatrixType.mass, nkel);
                        Ex = LocalMatrix(MatrixType.exotic2, nkel);
                        SR_GM_INGLOB(Velem, M, Direction.y);
                        for (i = 0; i < 3; i++)
                        {
                            for (j = 0; j < 3; j++)
                            {
                                t = GetValues(Velem[j], TypeOfSolution.Solution, Direction.y);
                                S_Y.F[Velem[i].globalNum] += Ex[i, j] * t;
                            }
                        }

                    }
                    S_Y.MSG();
                    flagY = 1;
                }
                return S_Y.Q;
            }
        }

        void SR_GM_INGLOB(point[] p, double[,] M, Direction d)   // добавление в глобальную матриуц (в зависимостри от напр. произв.)
        {
            int i, j, ind, k, s;
            if (d == Direction.x)
            {
                for (i = 0; i < 3; i++)
                    S_X.DI[p[i].globalNum] += M[i, i];

                for (k = 0; k < 3; k++)
                {
                    for (s = 0; s < 3; s++)
                    {
                        i = p[k].globalNum;
                        j = p[s].globalNum;
                        if (i < j)
                        {
                            for (ind = S_X.IG[j]; S_X.JG[ind] != i && ind <= S_X.IG[j + 1] - 1; ind++)
                                S_X.GGU[ind] += M[k, s];
                        }
                        else if (i > j)
                        {
                            for (ind = S_X.IG[i]; S_X.JG[ind] != j && ind <= S_X.IG[i + 1] - 1; ind++)
                                S_X.GGL[ind] += M[k, s];
                        }
                    }
                }
            }
            else
            {
                for (i = 0; i < 3; i++)
                    S_Y.DI[p[i].globalNum] += M[i, i];

                for (k = 0; k < 3; k++)
                {
                    for (s = 0; s < 3; s++)
                    {
                        i = p[k].globalNum;
                        j = p[s].globalNum;
                        if (i < j)
                        {
                            for (ind = S_Y.IG[j]; S_Y.JG[ind] != i && ind <= S_Y.IG[j + 1] - 1; ind++)
                                S_Y.GGU[ind] += M[k, s];
                        }
                        else if (i > j)
                        {
                            for (ind = S_Y.IG[i]; S_Y.JG[ind] != j && ind <= S_Y.IG[i + 1] - 1; ind++)
                                S_Y.GGL[ind] += M[k, s];
                        }
                    }
                }
            }
        }
    }

    public struct point
    {
        public double x;
        public double y;
        public int globalNum; // использовать только для вершин сетки/// глобальный номер узла в сетке
    }

    public enum ElemType { triangle, rectangle };//тип элемента:{треугольник, четырехугольник} 
    public enum BasisType { lin, quadr, cub }; // тип базиса:{линейный, квадратичный, тубический} 
    public enum MatrixType { mass, gest, exotic1, exotic2 };//тип матрицы:{массы, жесткости, особые} 
    public enum MaterialIdentifire { one, two, three };//номер подобласти/материала 
    public enum Direction { x, y }; //переменная, по которой дифференцируем
    public enum TypeOfSolution { Solution, Resultanta };//выдаем решение или результанту

    public class Element
    {
        public ElemType type;
        public BasisType btype;
        public point[] vertex; // массив вершин КЭ 
        public MaterialIdentifire material;

        public Element(ElemType t, MaterialIdentifire m, point[] v)
        {
            int i;
            type = t;
            material = m;
            if (type == ElemType.triangle)
            {
                vertex = new point[3];
                for (i = 0; i < 3; i++)
                {
                    vertex[i].x = v[i].x;
                    vertex[i].y = v[i].y;
                    vertex[i].globalNum = v[i].globalNum;
                }
            }
            else
            {
                vertex = new point[4];
                for (i = 0; i < 4; i++)
                {
                    vertex[i].x = v[i].x;
                    vertex[i].y = v[i].y;
                    vertex[i].globalNum = v[i].globalNum;
                }
            }
        }

    }


    public class Mesh : IEnumerable, IEnumerator
    {
        public Element[] elements;  //массив элементов 
        public double[] solution;   // массив весов 
        public Mesh(string FileNameMesh, string FileNameSolution)  //чтение сетки и решения из файлов  
        {
            FileStream f = null; 
            f = new FileStream(FileNameMesh, FileMode.Open, FileAccess.Read);
            StreamReader Reader = new StreamReader(f);
            string[] t;
            int i, N, N_elem, qwe, flag_of_vertex, j;
            string[] temp;
            ElemType type;
            point[] v; // массив вершин КЭ 
            MaterialIdentifire material;

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
            for (i = 0; i < 2; i++)
                t = Reader.ReadLine().Split(' ');
            t = Reader.ReadLine().Split(' ');
            N_elem = Convert.ToInt32(t[0]);
            elements = new Element[N_elem];
            for (i = 0; i < N_elem; i++)
            {
                temp = Reader.ReadLine().Split(' ');
                if (Convert.ToInt32(temp[1]) == 2)
                {
                    type = ElemType.triangle;
                    flag_of_vertex = 3;
                }
                else
                {
                    type = ElemType.rectangle;
                    flag_of_vertex = 4;
                }
                if (Convert.ToInt32(temp[3]) == 1)
                    material = MaterialIdentifire.one;
                else
                    if (Convert.ToInt32(temp[3]) == 2)
                        material = MaterialIdentifire.two;
                    else
                        material = MaterialIdentifire.three;
                v = new point[flag_of_vertex];
                for (j = 0; j < flag_of_vertex; j++)
                {
                    qwe = Convert.ToInt32(temp[6 + j]) - 1;
                    v[j].x = nodes[qwe].x;
                    v[j].y = nodes[qwe].y;
                    v[j].globalNum = qwe;
                }
                elements[i] = new Element(type, material, v);
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
        public double MaterialGamma(MaterialIdentifire type)
        {
            if (type == MaterialIdentifire.one)
                return 1.0;
            else
                if (type == MaterialIdentifire.two)
                    return 2.0;
                else
                    return 3.0;
        }
    }

    //**************************************************************
    // КЛАСС РЕШАТЕЛЬ
    //**************************************************************
    public class SLAE_MSG
    {
        public int Nmatrix, Nelem, Nig, Njg; // Nmatrix - кол-во точек, Nelem - кол-во элементов
        public int N, M;
        public int[] IG, JG;
        public double[] DI, GGU, GGL, F, F_LOC;
        public double[] Q;
        double[,] Dcoord_LOC2, M_LOC2;
        point[] Velem;

        public SLAE_MSG(string FileNameIG, string FileNameJG, int N) 
        {
            int i;
            string[] t;
            Nmatrix = N;
            Nelem = M;
            FileStream f = new FileStream(FileNameIG, FileMode.Open, FileAccess.Read);    //читаем IG
            StreamReader Reader = new StreamReader(f);
            Nig = Nmatrix + 1;
            IG = new int[Nig];
            for (i = 0; i < Nig; i++)
            {
                t = Reader.ReadLine().Split(' ');
                IG[i] = Convert.ToInt32(t[0]);
            }
            Reader.Close();

            f = new FileStream(FileNameJG, FileMode.Open, FileAccess.Read);     //читаем JG
            Reader = new StreamReader(f);
            Njg = IG[Nmatrix];
            JG = new int[Njg];
            for (i = 0; i < Njg; i++)
            {
                t = Reader.ReadLine().Split(' ');
                JG[i] = Convert.ToInt32(t[0]);
            }
            Reader.Close();

            DI = new double[Nmatrix];
            GGU = new double[Njg];
            GGL = new double[Njg];
            F = new double[Nmatrix];
            Q = new double[Nmatrix];
            F_LOC = new double[3];
            M_LOC2 = new double[3, 3];   //!!!
            Dcoord_LOC2 = new double[3, 3];

            Velem = new point[3];
        }

        void mult(double a, double[] x, double[] y)
        {
            for (int i = 0; i < Nmatrix; i++)
                y[i] = x[i] + a * y[i];
        }

        void summult(double a, double[] x, double[] y)
        {
            for (int i = 0; i < Nmatrix; i++)
                y[i] += a * x[i];
        }

        void copy(double[] x, double[] y)
        {
            for (int i = 0; i < Nmatrix; i++)
                y[i] = x[i];
        }

        double norma(double[] x)
       {
           double res = 0.0;
            for (int i=0; i<Nmatrix; i++)
             res += Math.Pow(x[i],2);
           res = Math.Sqrt(res);
           return res;
       }

        double scal(double[] x, double[] y)
       {
           double res = 0.0;
            for (int i=0; i<Nmatrix; i++)
             res += x[i] * y[i];
           return res;
       }

        void mult_A(double[] x, double[] y)
        {
            int j, ig0, ig1;

            for (int i = 0; i < Nmatrix; i++)
            {
                ig0 = IG[i];
                ig1 = IG[i + 1];
                y[i] = DI[i] * x[i];
                for (int k = ig0; k < ig1; k++)
                {
                    j = JG[k];
                    y[i] += GGL[k] * x[j];
                    y[j] += GGU[k] * x[i];
                }
            }
        }

        void mult_f_A(double[] x, double[] y)
        {
            int j, i, ig0, ig1;

            for (i = 0; i < Nmatrix; i++)
            {
                ig0 = IG[i];
                ig1 = IG[i + 1];
                y[i] = F[i] - DI[i] * x[i];
                for (int k = ig0; k < ig1; k++)
                {
                    j = JG[k];
                    y[i] -= GGL[k] * x[j];
                    y[j] -= GGU[k] * x[i];
                }
            }
        }

        public void MSG()
        {
            int i;
             int maxit=1, MAX=10000;
             double A_RESH, B_RESH, T1, nev = 1.0, f, e=1e-15;
             double[] R_RESH, Z_RESH, T2, HL;
             
             R_RESH = new double [Nmatrix];  //#######
             Z_RESH = new double [Nmatrix];  //#######
             T2 = new double[Nmatrix];  //#######
             HL = new double[Nmatrix];  //#######

          for(i=0; i<Nmatrix-1; i++)
           Q[i]=1.0;
          Q[Nmatrix-1]=1.0;
     
         mult_f_A(Q,R_RESH);   // f-Ax0  //#######
         copy (R_RESH,Z_RESH);
         f = norma(F);                // ||pr||  //#######

          while(nev>e && maxit <= MAX)  //#######
         {
          T1 = scal(R_RESH,R_RESH);                  // (rk-1, rk-1) #######
          mult_A(Z_RESH,T2);   // Azk-1  //#######
          A_RESH = T1/scal(T2,Z_RESH);               // Ak #######
          summult (A_RESH,Z_RESH,Q);                 // xk #######
          summult (-A_RESH,T2,R_RESH);               // rk #######
          B_RESH = scal(R_RESH,R_RESH)/T1;                  // Bk #######
          mult (B_RESH,R_RESH,Z_RESH);                    // zk #######
          nev = norma(R_RESH)/f;                // otnos nev  //#######
          maxit++;  //#######
         }
        }
    }
}