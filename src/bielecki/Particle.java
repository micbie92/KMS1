package bielecki;

import java.util.Random;

public class Particle {

    public static double kb = 0.00831;
    
    
    private double Enumber=0;
    private Point E;
    private Point r;
    private Point p;

    public Particle(Point r, double T) {
        this.r = r;
        
        Random rand = new Random();

        double lambda = rand.nextDouble();
        double Ex = (-1 / 2.) * kb * T * Math.log(lambda);
        double px = Math.sqrt(2 * Argon.m * Ex);
        if (rand.nextDouble() < 0.5)
            px = -px;

        lambda = rand.nextDouble();
        double Ey = (-1 / 2.) * kb * T * Math.log(lambda);
        double py = Math.sqrt(2 * Argon.m * Ey);
        if (rand.nextDouble() < 0.5)
            py = -py;

        lambda = rand.nextDouble();
        double Ez = (-1 / 2.) * kb * T * Math.log(lambda);
        double pz = Math.sqrt(2 * Argon.m * Ez);
        if (rand.nextDouble() < 0.5)
            pz = -pz;

        Enumber = Ex + Ey + Ez;
        E = new Point(Ex, Ey, Ez);
        p = new Point(px, py, pz);
    }

    public String toStringg() {
        return String.format("%7f\t%7f\t%7f\n",r.x,r.y,r.z);
    }

    public Point getE() {
        return E;
    }

    public double getEnumber() {
        return Enumber;
    }

    public Point getR() {
        return r;
    }

    public void setR(Point r) {
        this.r = r;
    }

    public Point getP() {
        return p;
    }

    public void setP(Point p) {
        this.p = p;
        E.x = p.x*p.x/(2*Argon.m);
        E.y = p.y*p.y/(2*Argon.m);
        E.z = p.z*p.z/(2*Argon.m);
        Enumber = E.x + E.y + E.z;
        
    }

}