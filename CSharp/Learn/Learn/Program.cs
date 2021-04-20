using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Learn
{
	class Rectangle
	{
		double length;
		double width;

		public Rectangle(double length, double width)
		{
			this.length = length;
			this.width = width;
		}

		public void SetValue(double length, double width)
		{
			this.length = length;
			this.width = width;
		}

		public double GetLength()
		{
			return length;
		}

		public double GetWidth()
		{
			return width;
		}

		public void PrintValue()
		{
			Console.WriteLine("Length is {0}, Width is {1}.", length, width);
		}
	}

	class Program
	{
		static void Main(string[] args)
		{
			int[,,] mat = new int[4, 3, 2];

			for (int i = 0; i != 4; ++i)
				for (int j = 0; j != 3; ++j)
					for (int k = 0; k != 2; ++k)
					mat[i, j, k] = i + j + k + 4;

			foreach (var i in mat)
				Console.WriteLine(i);

			Console.ReadKey();
		}
	}
}
