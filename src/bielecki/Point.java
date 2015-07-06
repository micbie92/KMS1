package bielecki;

public class Point {

    double x;
    double y;
    double z;

    Point() {
        x = 0;
        y = 0;
        z = 0;
    }
    
    Point(double X, double Y, double Z) {
        x = X;
        y = Y;
        z = Z;
    }

    
    public Point add(Point p) {
        return new Point(this.x + p.x, this.y + p.y, this.z + p.z);
    }
    
    public Point multiply(double c) {
        return new Point(this.x * c, this.y * c, this.z * c);
    }

    public static double distanceTwo(Point p1, Point p2) {
        double dsquare = (p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y) + (p1.z - p2.z) * (p1.z - p2.z);
        return Math.sqrt(dsquare);
    }
  
    public static double distance(Point p) {
        double lengthSquare = p.x * p.x + p.y * p.y + p.z * p.z;
        return Math.sqrt(lengthSquare);
    }
    
    static Point minus(Point p1, Point p2) {
        return new Point(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z);
    }    

}