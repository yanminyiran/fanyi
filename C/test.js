//alert("hello")
/*private static var gg = 9.8 * m0 / S * 0.001;//导线荷载，
        private static var alpha = 113.68;//导线应力系数，固定参数

        static var D = 15.07;//样例导线直径,mm
        static var S = 134.49;//样例导线横截面面积,mm2
        static var m0 = 466.8;//样例导线单位长度质量,kg/km
        static var KWind = 1.1;//风载体型系数
        static var Awind = 1;//风压不均匀系数
        static var tempr = 20;//初始温度*/
var D = 15.07;
var S=134.49
var Awind = 1;
var tempr = 20;
var KWind = 1.1;
var alpha = 0.00002;
var m0 = 466.8;

var E = 78400;


function Point(x,y,z,sort){
  this.x =x;
  this.y =y;
  this.z =z;
  this.sort=sort;
}

 //输出 AB间距，此函数用于FinalCal
 function Distance(p1, p2,zAvailable)
        {
           if (zAvailable)
            {
                return Math.sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y) + (p1.z - p2.z) * (p1.z - p2.z));
            }
            else
            {
                return Math.sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
            }
        }
		//测试
//输入 点p1、p2，是否考虑Z方向距离zAvailable
/**
var p1=new Point(0,0,0);
 var p2=new Point(3,4,4);
alert(Distance(p1,p2,0));
*/

function getAngle(x1, y1,  x2,  y2, windD)//计算风向与两杆塔间线路方向的夹角
{
    var angle = 0.0;
    var va_x = x2 - x1;
    var va_y = y2 - y1;
    var vb_x = Math.cos(windD / 180 * Math.PI);
    var vb_y = Math.sin(windD / 180 * Math.PI);
    var productValue = (va_x * vb_x) + (va_y * vb_y);  
    var va_val = Math.sqrt(va_x * va_x + va_y * va_y);       
    var vb_val = Math.sqrt(vb_x * vb_x + vb_y * vb_y);        
    var cosValue = productValue / (va_val * vb_val);
    if (cosValue < -1 && cosValue > -2)
    { cosValue = -1; }
    else if (cosValue > 1 && cosValue < 2)
    { cosValue = 1; }
    angle = Math.acos(cosValue) * 180 / Math.PI;
    return angle;
	}
//测试
alert(getAngle(4, 4,5,5,0));

function TripleX(a,b,c,d){
	var result = 0.0;
try
   {
	var r1 = 0.5 * d / a - b * c / a / a / 6 + b * b * b / a / a / a / 27;
	var r2 = c / a / 3 - b * b / a / a / 9;
	result = result - b / 3 / a;
	var rmid = Math.sqrt(r1 * r1 + r2 * r2 * r2);
	if (-r1 + rmid < 0 && -r1 - rmid < 0)
    {
        result = result - Math.pow(Math.abs(-r1 + rmid), 1.0*d / (3.0*d)) - Math.pow(Math.abs(-r1 - rmid), 1.0*d / 3.0*d);
        return result;
    }
    else if (-r1 + rmid > 0 && -r1 - rmid < 0)
    {
        result = result + Math.pow(Math.abs(-r1 + rmid), 1.0*d / 3.0*d) - Math.Pow(Math.Abs(-r1 - rmid), 1.0*d / 3.0*d);
        return result;
    }
    else if (-r1 + rmid < 0 && -r1 - rmid > 0)
    {
        result = result - Math.pow(Math.abs(-r1 + rmid), 1.0*d / 3.0*d) + Math.Pow(Math.abs(-r1 - rmid), 1.0*d / 3.0*d);
        return result;
    }
    else if (-r1 + rmid > 0 && -r1 - rmid > 0)
    {
        result = result + Math.pow(Math.abs(-r1 + rmid), 1.0*d / 3.0*d) + Math.pow(Math.abs(-r1 - rmid), 1.0*d / 3.0*d);
        return result;
    }
	}
            catch (err)
            {
                return 0.0;
            }
	return result;
}




function JudgeDirect(x1, y1,  x2, y2, windD)
{
    var minus = 1;
    var r = Math.sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
    var x0 = (x2 + x1) / 2 + r * Math.cos(windD / 180 * Math.PI);
    var y0 = (y1 + y2) / 2 + r * Math.sin(windD / 180 * Math.PI);
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




function sortNumber(a,b)
{
	return b-a
}

function LowCal(origin,p1, p2)
        {
            var affected = new Array();
            for (var i = 0; i < origin.length; i++)
            {
                var temp = new Point();
                temp.x = origin[i].x;
                temp.y = origin[i].y;
                temp.z = p1.z + (origin[i].x - p1.x) / (p2.x - p1.x) * (p2.z - p1.z) - origin[i].z;
                affected[i]=temp;
            }
            lowPoint = affected.sort(sortNumber);
            return lowPoint[0];
        }




//输入 导线点集origin，起始点T1，终止点T2，目标温度，数据生成时温度
        //输出 生成时温度应力与目标温度应力的比值result，此函数用于FinalCal
       function TempCal(origin, t1, t2, tempSet, tempNow)
        {
            var result = 0.0;
            lowPoint = LowCal(origin, t1, t2);
            var gObserve = 9.8 * m0 / S * 0.001 * Distance(t1, lowPoint, false) * Distance(t2, lowPoint, false) / 2 / lowPoint.z;
            var b = gObserve - E * 9.8 * m0 / S * 0.001 * 9.8 * m0 / S * 0.001 * Distance(t1, lowPoint, false) * Distance(t2, lowPoint, false) / 6 / gObserve / gObserve - alpha * E * (tempNow - tempSet);
            var d = E * 9.8 * m0 / S * 0.001 * 9.8 * m0 / S * 0.001 * Distance(t1, lowPoint, false) * Distance(t2, lowPoint, false) / 6;
            var gTemp = TripleX(1, b, 0, d);
            result = gObserve / gTemp;
            return result;
        }







//输入 导线点集origin，起始点T1，终止点T2，风向，风速，冰厚，目标温度，数据生成时温度
        //输出 最终计算得到的点集
        function FinalCal(origin, t1,  t2,  windD, windV, ice,  tempr, tempNow)
        {
            var g1 = 9.8 * m0 / S * 0.001;
            var g2 = 27.708 * (ice + D) * ice / S * 0.001;
            var angleWind = getAngle(t1.x, t1.y, t2.x, t2.y, windD);
            windV = windV * Math.sin(angleWind / 180 * Math.PI);
            var g5 = 0.6125 * KWind * Awind * (2 * ice + D) * windV * windV / S * 0.001;
            var g = Math.sqrt((g1 + g2) * (g1 + g2) + g5 * g5);
            var sita = Math.atan(g5 / (g1 + g2)) * 180 / Math.PI;//风偏转角
            var sigma = Math.atan((t2.y - t1.y) / (t2.x - t1.x)) * 180 / Math.PI;//线路走向与坐标轴偏角
            var dire = JudgeDirect(t1.x, t1.y, t2.x, t2.y, windD);

            var Temprindex = new Array();
            for (var i = -20; i < 51; i++)
            {
                Temprindex[i]=TempCal(origin, t1, t2, i, 50);
            }
            var tempNum = Temprindex[tempr + 20] / Temprindex[tempNow + 20];
            var affected = new Array();
            for (var i = 0; i < origin.length; i++)
            {
                affected[i]=((t1.z + (origin[i].x - t1.x) / (t2.x - t1.x) * (t2.z - t1.z) - origin[i].z) * g / g1) * tempNum;
            }
            var result = new Array();
            for (var i = 0; i < origin.length; i++)
            {
                var temp = new Point();
                temp.sort = 16;
                temp.x = origin[i].x + affected[i] * Math.sin(sita / 180 * Math.PI) * Math.sin(sigma / 180 * Math.PI) * dire;
                temp.y = origin[i].y + affected[i] * Math.sin(sita / 180 * Math.PI) * Math.cos(sigma / 180 * Math.PI) * dire;
                temp.z = origin[i].z - affected[i] * (1 - Math.cos(sita / 180 * Math.PI));
                //temp.x = affected[i] * Math.Cos(sigma / 180 * Math.PI) * Math.Sin(sita / 180 * Math.PI) ;
                //temp.y = affected[i] * Math.Sin(sigma / 180 * Math.PI) * Math.Sin(sita / 180 * Math.PI) ;
                //temp.z = affected[i] * (1 - Math.Cos(sita / 180 * Math.PI));
                result[i]=temp;
            }           
            return result;
        }



var origin=new Array();
origin[0]=new Point(3,1,2);
origin[1]=new Point(2,2,3);
origin[2]=new Point(0,1,4);
origin[3]=new Point(3,2,5);
var p1=new Point(0,0,0);
var p2=new Point(3,4,4);

var finalResult=FinalCal(origin, p1,  p2,  3, 5, 4,  20, 25);

for( var i=0;i<finalResult.length;i++){
document.write("point"+i+".x="+finalResult[i].x+'\t');
document.write("point"+i+".y="+finalResult[i].y+'\t');
document.write("point"+i+".z="+finalResult[i].z+"<br />");
}

