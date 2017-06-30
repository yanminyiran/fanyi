using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace DrawThrowLine4
{
    class WeatherControlClass
    {
        static double D = 15.07;//样例导线直径,mm
        static double S = 134.49;//样例导线横截面面积,mm2
        static double m0 = 466.8;//样例导线单位长度质量,kg/km
        static double KWind = 1.1;//风载体型系数
        static double Awind = 1;//风压不均匀系数
        static double tempr = 20;//初始温度
        static double E = 78400;//弹性模量,N/mm2
        static double aerfa = 0.00002;//金属线膨胀系数,1/°C
        static double alpha = 113.68;//导线应力系数,N/mm2

        public class Point
        {
            public double x { get; set; }
            public double y { get; set; }
            public double z { get; set; }
            public int sort { get; set; }
        }

        //输入 点A、B，是否考虑Z方向距离zAvailable
        //输出 AB间距，此函数用于FinalCal
        private static double Distance(Point a, Point b, bool zAvailable)
        {
            if (zAvailable)
            {
                return Math.Sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z));
            }
            else
            {
                return Math.Sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
            }
        }

        //输入 两点1、2的x、y方向坐标，风向windD
        //输出 风向与两杆塔间线路方向的夹角angle，此函数用于FinalCal
        public static double getAngle(double x1, double y1, double x2, double y2, int windD)
        {
            double angle = 0.0;
            double va_x = x2 - x1;
            double va_y = y2 - y1;
            double vb_x = Math.Cos((double)windD / 180 * Math.PI);
            double vb_y = Math.Sin((double)windD / 180 * Math.PI);
            double productValue = (va_x * vb_x) + (va_y * vb_y);
            double va_val = Math.Sqrt(va_x * va_x + va_y * va_y);
            double vb_val = Math.Sqrt(vb_x * vb_x + vb_y * vb_y);
            double cosValue = productValue / (va_val * vb_val);
            if (cosValue < -1 && cosValue > -2)
            { cosValue = -1; }
            else if (cosValue > 1 && cosValue < 2)
            { cosValue = 1; }
            angle = Math.Acos(cosValue) * 180 / Math.PI;
            return angle;
        }

        //输入 ax3+bx2+cx+d=0方程的系数a、b、c、d
        //输出 方程实根result，此函数用于TempCal
        public static double TripleX(double a, double b, double c, double d)
        {
            double result = 0.0;
            try
            {
                double r1 = 0.5 * d / a - b * c / a / a / 6 + b * b * b / a / a / a / 27;
                double r2 = c / a / 3 - b * b / a / a / 9;
                result = result - b / 3 / a;
                double rmid = Math.Sqrt(r1 * r1 + r2 * r2 * r2);
                if (-r1 + rmid < 0 && -r1 - rmid < 0)
                {
                    result = result - Math.Pow(Math.Abs(-r1 + rmid), 1.0d / 3.0d) - Math.Pow(Math.Abs(-r1 - rmid), 1.0d / 3.0d);
                    return result;
                }
                else if (-r1 + rmid > 0 && -r1 - rmid < 0)
                {
                    result = result + Math.Pow(Math.Abs(-r1 + rmid), 1.0d / 3.0d) - Math.Pow(Math.Abs(-r1 - rmid), 1.0d / 3.0d);
                    return result;
                }
                else if (-r1 + rmid < 0 && -r1 - rmid > 0)
                {
                    result = result - Math.Pow(Math.Abs(-r1 + rmid), 1.0d / 3.0d) + Math.Pow(Math.Abs(-r1 - rmid), 1.0d / 3.0d);
                    return result;
                }
                else if (-r1 + rmid > 0 && -r1 - rmid > 0)
                {
                    result = result + Math.Pow(Math.Abs(-r1 + rmid), 1.0d / 3.0d) + Math.Pow(Math.Abs(-r1 - rmid), 1.0d / 3.0d);
                    return result;
                }
            }
            catch (Exception ex)
            {
                return 0.0;
            }
            return result;
        }

        //输入 两点1、2的x、y方向坐标，风向windD
        //输出 风导致两点连线的偏移方向（正负，以minus表示），此函数用于FinalCal
        public static int JudgeDirect(double x1, double y1, double x2, double y2, int windD)
        {
            int minus = 1;
            double r = Math.Sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
            double x0 = (x2 + x1) / 2 + r * Math.Cos(windD / 180 * Math.PI);
            double y0 = (y1 + y2) / 2 + r * Math.Sin(windD / 180 * Math.PI);
            if (x0 > x2 && x0 > x1)
            {
                minus = -1;
            }
            else if (y0 > y2 && y0 > y1)
            {
                minus = -1;
            }
            return minus;
        }

        //输入 导线点集origin，起始点T1，终止点T2
        //输出 导线最大弧垂点lowPoint，此函数用于TempCal
        private static Point LowCal(List<Point> origin, Point t1, Point t2)
        {
            List<Point> affected = new List<Point>();
            for (int i = 0; i < origin.Count; i++)
            {
                Point temp = new Point();
                temp.x = origin[i].x;
                temp.y = origin[i].y;
                temp.z = t1.z + (origin[i].x - t1.x) / (t2.x - t1.x) * (t2.z - t1.z) - origin[i].z;
                affected.Add(temp);
            }
            Point lowPoint = affected.OrderByDescending(p => p.z).First();
            return lowPoint;
        }

        //输入 导线点集origin，起始点T1，终止点T2，目标温度，数据生成时温度
        //输出 生成时温度应力与目标温度应力的比值result，此函数用于FinalCal
        public static double TempCal(List<Point> origin, Point t1, Point t2, double tempSet, double tempNow)
        {
            double result = 0.0;
            Point lowPoint = LowCal(origin, t1, t2);
            double gObserve = 9.8 * m0 / S * 0.001 * Distance(t1, lowPoint, false) * Distance(t2, lowPoint, false) / 2 / lowPoint.z;
            double b = gObserve - E * 9.8 * m0 / S * 0.001 * 9.8 * m0 / S * 0.001 * Distance(t1, lowPoint, false) * Distance(t2, lowPoint, false) / 6 / gObserve / gObserve - aerfa * E * (tempNow - tempSet);
            double d = E * 9.8 * m0 / S * 0.001 * 9.8 * m0 / S * 0.001 * Distance(t1, lowPoint, false) * Distance(t2, lowPoint, false) / 6;
            double gTemp = TripleX(1, b, 0, d);
            result = gObserve / gTemp;
            return result;
        }

        //输入 导线点集origin，起始点T1，终止点T2，目标温度，数据生成时温度，目标状态下比载，标准比载
        //输出 生成时温度应力与目标温度应力的比值result，此函数用于FinalCal
        public static double TempCal2(List<Point> origin, Point t1, Point t2, double tempSet, double tempNow, double GSet, double GNow)
        {
            double result = 0.0;
            double b = alpha - E * GSet * GSet * Distance(t1, t2, false) * Distance(t1, t2, false) / 24 / alpha / alpha - aerfa * E * (tempNow - tempSet);
            double d = E * GNow * GNow * Distance(t1, t2, false) * Distance(t1, t2, false) / 24;
            result = alpha / TripleX(1, -b, 0, -d);
            return result;
        }

        //输入 导线点集origin，起始点T1，终止点T2，风向，风速，冰厚，目标温度，数据生成时温度
        //输出 最终计算得到的点集
        public static List<Point> FinalCal(List<Point> origin, Point t1, Point t2, int windD, double windV, int ice, int tempr, int tempNow)
        {
            double g1 = 9.8 * m0 / S * 0.001;
            double g2 = 27.708 * (ice + D) * ice / S * 0.001;
            double angleWind = getAngle(t1.x, t1.y, t2.x, t2.y, windD);
            windV = windV * Math.Sin(angleWind / 180 * Math.PI);
            double g5 = 0.6125 * KWind * Awind * (2 * ice + D) * windV * windV / S * 0.001;
            double g = Math.Sqrt((g1 + g2) * (g1 + g2) + g5 * g5);
            double sita = Math.Atan(g5 / (g1 + g2)) * 180 / Math.PI;//风偏转角
            double sigma = Math.Atan((t2.y - t1.y) / (t2.x - t1.x)) * 180 / Math.PI;//线路走向与坐标轴偏角
            int dire = JudgeDirect(t1.x, t1.y, t2.x, t2.y, windD);
            //List<double> Temprindex = new List<double>();
            //for (int i = -20; i < 51; i++)
            //{
            //    Temprindex.Add(WeatherControlClass.TempCal(origin, t1, t2, i, 50));
            //}
            //double tempNum = Temprindex[tempr + 20] / Temprindex[tempNow + 20];
            double alphaNum = WeatherControlClass.TempCal2(origin, t1, t2, tempr, tempNow, g1, g) * g / g1;
            List<double> affected = new List<double>();
            for (int i = 0; i < origin.Count; i++)
            {
                affected.Add((t1.z + (origin[i].x - t1.x) / (t2.x - t1.x) * (t2.z - t1.z) - origin[i].z) * alphaNum);
            }
            List<Point> result = new List<Point>();
            for (int i = 0; i < origin.Count; i++)
            {
                Point temp = new Point();
                temp.sort = 16;
                temp.x = origin[i].x + affected[i] * Math.Sin(sita / 180 * Math.PI) * Math.Sin(sigma / 180 * Math.PI) * dire;
                temp.y = origin[i].y + affected[i] * Math.Sin(sita / 180 * Math.PI) * Math.Cos(sigma / 180 * Math.PI) * dire;
                temp.z = origin[i].z - affected[i] * Math.Cos(sita / 180 * Math.PI);
                //temp.x = affected[i] * Math.Cos(sigma / 180 * Math.PI) * Math.Sin(sita / 180 * Math.PI) ;
                //temp.y = affected[i] * Math.Sin(sigma / 180 * Math.PI) * Math.Sin(sita / 180 * Math.PI) ;
                //temp.z = affected[i] * (1 - Math.Cos(sita / 180 * Math.PI));
                result.Add(temp);
            }
            //result.OrderByDescending(p => p.x);
            return result;
        }
    }
}
