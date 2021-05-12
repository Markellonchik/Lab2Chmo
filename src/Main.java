public class Main {
	static int m = 10;
	static int n = 4;
	public static void main(String[] args) {
		solve2();
	}
	
	static void solve2() {
		System.out.println("Task 2:");
		
		double x[] = new double[m];
		double y[] = new double[m];
		double last = 0;
		for(int i = 0; i < m; ++i) {
			x[i] = last;
			last += 0.2;
		}
		
		System.out.print("x: [");
		for(double xi : x) {
			System.out.print(xi + ",");
		}
		System.out.println("]");
		
		for(int i = 0; i < m; ++i) {
			y[i] = f(x[i]);
		}
		
		x = new double[]{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1};
		y = new double[]{0.29, 0.17, 0.05, -0.07, -0.19, -0.32, -0.5, -0.71, -0.8, -1};
		
		printArray(x);
		printArray(y);
		
		double[][] G = new double[n + 1][n + 1];
		for(int i = 0; i <= n; ++i) {
			for(int j = 0; j <= n; ++j) {
				G[i][j] = 0;
				for(int k = 0; k < m; ++ k) {
					G[i][j] += Math.pow(x[k], i) * Math.pow(x[k], j);
				}
			}
		}
		
		double fG[] = new double[n + 1];
		for(int i = 0; i <= n; ++i) {
			fG[i] = 0;
			for(int k = 0; k < m; ++k) {
				fG[i] += (i == 0 ? 1 : Math.pow(x[k], i)) * y[k];
			}
		}
		
		System.out.println("Матрица Грама");
		printMatrix(G);
		
		printArray(fG);
		
		elimination(G, fG);
		double[] q = backSubstitution(G, fG);
		
		System.out.println("Коэфициенты арифметического приближения");
		printArray(q);
		
		double dif = vidhilennya(q, x, y);
		System.out.println("Vidhilennya: " + dif);
	}
	
	public static double okruglenie(double x) {
		int n = (int) Math.round(x);
		return (double) n / 100;
	}
	
	public static double vidhilennya(double[] q, double[] x, double[] f) {
		double res = 0;
		for(int k = 0; k < m; ++k) {
			double r = f[k];
			for(int i = 0; i <= n; ++i) {
				r -= (i == 0 ? 1 : Math.pow(x[k], i)) * q[i];
			}
			res += r * r;
		}
		return res;
	}

	public static void printMatrix(double[][] A) {
		System.out.print("[");
		for(int i = 0; i < A.length; ++i) {
			for(int j = 0; j < A[i].length; ++j) {
				if(j > 0) System.out.print(", ");
				System.out.print(A[i][j]);
			}
			if(i + 1 == A.length) System.out.print("]");
			System.out.println(";");
		}
	}
	public static void printArray(double[] a) {
		System.out.print("[");
		for(int i = 0; i < a.length; ++i) {
			if(i > 0) System.out.print(", ");
			System.out.print(a[i]);
		}
		System.out.println("]");
	}
	
	public static void elimination(double[][] A, double[] b) {
        for(int k = 0; k <= n; k++) {
            for(int i = k+1; i <= n; i++) {
                double l = A[i][k] / A[k][k];
                A[i][k] = 0;
                
                for(int j = k+1; j <= n; j++) {
                    A[i][j] = A[i][j] - l * A[k][j];
                }
                b[i] = b[i] - l * b[k];
            }
        }
    }
	
	public static double[] backSubstitution(double[][] A, double[] b) {
		double x[] = new double[n + 1];
		for(int i = 0; i <= n; ++i) x[i] = 0;
		
        x[n] = b[n] / A[n][n];
        for(int i = n - 1; i >= 0; i--) {
            x[i] = (b[i] - solve(i, A, x)) / A[i][i];
        }
        return x;
    }
	
	public static double solve(int i, double A[][], double[] x) {
        double result = 0.0;
        for(int j = i; j <= n; j++)
            result += A[i][j] * x[j];
        return result;
    }
	
	static double f(double x) {
		return Math.sin(x) - Math.log(x);
	}
}
