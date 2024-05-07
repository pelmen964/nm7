using System;
using System.IO;
using AngouriMath;
using AngouriMath;
using AngouriMath.Extensions;
using static System.Console;
using static AngouriMath.MathS;

namespace nm7
{
    internal class Program
    {
        public static void Main(string[] args)
        {
            Task task = new Task("TaskInput.txt");
            Out.WriteLine(task.IntegrateWithLefthanded());
            Out.WriteLine(task.IntegrateWithRighthanded());
            Out.WriteLine(task.IntegrateWithTrapeze());
            // double[] a = new[] { 0, 0.1, 0.3, 0.5, 0.8, 1.2, 2, 3 };
            // foreach (var j1 in a)
            // {
            //     Out.Write(task.F(j1)+" ");
            // }
            
        }
    }
}