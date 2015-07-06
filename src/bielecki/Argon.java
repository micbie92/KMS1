package bielecki;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

public class Argon {

    static double a = 0;
    static int n = 0;
    static int aaa = 0;
    static double T = 0;
    static double tau = 0;
    static double L = 0;
    static double R = 0;
    static double m = 0;
    static double eps = 0;
    static double f = 0;
    
    static int So = 0;
    static int Sd = 0;
    static int Sout = 0;
    static int Sxyz = 0;

    static Point b1;
    static Point b2;
    static Point b3;

    static int N;
    static double P;
    static double Ek;
    static double H;
    
    static double time=0;
    static double Tsum=0;
    static double Esum=0;
    static double Psum=0;
    static double Hsum=0;

    static List<Particle> particleList;
    static double[][] potentialParticles;
    static double[] potentialSp;
    static double potentialSum;

    static Point[] forceParticles;
    static Point[] forceSp;
    static Point[] forceSum;

    public static void main(String[] args) throws IOException {

        loadData("dane");
        initVectors();
        createParticles();
        
        potentialParticles = new double[N][N];
        potentialSp = new double[N];
        calculatePotentialSum();
       /* 
        File fileOutt = new File("pot");
        FileOutputStream fopp = new FileOutputStream(fileOutt);
        
        for(int i=0;i<800;i++)
        {
        	a +=0.001;
        	initVectors();
            createParticles();
            calculatePotentialSum();
            fopp.write(String.format("%f\t%f\n", a,potentialSum).getBytes());
        	
        }
        fopp.close();
*/
        forceParticles = new Point[N];
        forceSp = new Point[N];
        forceSum = new Point[N];   
        calculateForceSum();
        calculatePreasure();
        calculateTemperature();
        calculateTotalEnergy();
        checkMomentum();
        checkAverageEnergy();
        normalizeMomentum();
        checkAverageEnergy();
        checkMomentum(); 
        for (int s=1;s<So;s++){
        	calculateP();
        	calculateR();
        	calculateForceSum();
        	calculatePotentialSum();
            calculatePreasure();
            calculateP();
        	checkAverageEnergy();
        	checkMomentum();
            calculateTemperature();
            calculateTotalEnergy();
        }
        File fileOutt = new File("wyniki");
        FileOutputStream fopp = new FileOutputStream(fileOutt);
        saveData(fopp);
        
        File fileOut = new File("XYZ.dat");
        FileOutputStream fop = new FileOutputStream(fileOut);
        saveXYZ(fop);
        
        
      //  File fileOutt = new File("energy");
      //  FileOutputStream fopp = new FileOutputStream(fileOutt);
        int n=1;
        int m=1;
        for (int s=1;s<Sd;s++){
        	calculateP();
        	calculateR();
        	calculateForceSum();
        	calculatePotentialSum();
            calculatePreasure();
            calculateP();
        	checkAverageEnergy();
        	checkMomentum();
            calculateTemperature();
            calculateTotalEnergy();
            time+=tau;
            Tsum+=T;
            Esum+=Ek;
            Psum+=P;
            Hsum+=H;
        	if(s==n*Sxyz){
        		saveXYZ(fop);
        		n++;
        	}
        	if (s==m*Sout){
        		saveData(fopp);
        		//fopp.write(String.format("%f\t%f\n", time,H).getBytes());
        		m++;
        	}
        }
        
       // fop.close();
        fopp.close();
        
        double Tsr=Tsum/(Sd-1);
        double Esr=Esum/(Sd-1);
        double Psr=Psum/(Sd-1);
        double Hsr=Hsum/(Sd-1);
        System.out.println(String.format("Tsr:%f  ;Eksr:%f  ;Psr:%f  ;Ecsr:%f",Tsr,Esr,Psr,Hsr));
       
        System.out.println("Zakończono wykonywanie programu");

    }
    
    
    
// FUNKCJE
    
    private static void loadData(String fileName) {
        Scanner scan;
        File fileIn = new File(fileName);
        try {
            scan = new Scanner(fileIn);
            a = scan.nextDouble();
            n = scan.nextInt();
            T = scan.nextDouble();
            tau = scan.nextDouble();
            L = scan.nextDouble();
            R = scan.nextDouble();
            m = scan.nextDouble();
            eps = scan.nextDouble();
            f = scan.nextDouble();
            So = scan.nextInt();
            Sd = scan.nextInt();
            Sout = scan.nextInt();
            Sxyz = scan.nextInt();
            
        } catch (Exception e1) {
            System.out.println("Invalid input file");
        }

        System.out.println(String.format("Wczytane dane:\na=%.3f, n=%d, T=%f, tau=%.3f, L=%.2f, R=%.3f, m =%.1f, eps=%.1f, f=%.1f So=%d, Sd=%d, Sout=%d, Sxyz=%d", a, n, T, tau, L, R, m, eps, f, So, Sd, Sout, Sxyz));
    }

    private static void saveData(FileOutputStream fopp) {
        try {
             //fopp.write(String.format("t:%.2f; Ec:%.3f; Ek:%.3f; Ep:%.3f; Temp:%.1f; P:%.1f\n",time,H,Ek,potentialSum,T,P).getBytes());
             fopp.write(String.format("%f\t%f\t%f\t%f\t%f\t%f\n",time,H,Ek,potentialSum,T,P).getBytes());
        } catch (IOException e) {
            System.out.println("File exception");
        }
    }
    
    private static void saveXYZ(FileOutputStream fop) {
        try {
        	aaa++;
        	System.out.println(String.format("XYZ %d", aaa));
            for (Particle p : particleList) {
                fop.write(p.toStringg().getBytes());
            }
            fop.write("\n\n".getBytes());
        } catch (IOException e) {
            System.out.println("File exception");
        }
    }
    

    private static void initVectors() {
        b1 = new Point(a, 0, 0);
        b2 = new Point(a / 2, a * Math.sqrt(3) / 2, 0);
        b3 = new Point(a / 2, a * Math.sqrt(3) / 6, a * Math.sqrt(2. / 3));
    }

    private static void createParticles() {
        particleList = new ArrayList<Particle>();
        double c = (n - 1) / 2.;

        for (int i1 = 0; i1 < n; i1++) {
            for (int i2 = 0; i2 < n; i2++) {
                for (int i3 = 0; i3 < n; i3++) {
                    double x = (i1 - c) * b1.x + (i2 - c) * b2.x + (i3 - c) * b3.x;
                    double y = (i1 - c) * b1.y + (i2 - c) * b2.y + (i3 - c) * b3.y;
                    double z = (i1 - c) * b1.z + (i2 - c) * b2.z + (i3 - c) * b3.z;
                    Point r = new Point(x, y, z);
                    particleList.add(new Particle(r, T));
                }
            }
        }
        N = particleList.size();
    }

    private static void calculatePotentialSp() {
        for (int i = 0; i < N; i++) {
            Particle p = particleList.get(i);
            double ri = Point.distance(p.getR());
            if (ri < L) {
                potentialSp[i] = 0;
            } else {
                potentialSp[i] = 0.5 * f * (ri - L) * (ri - L);
            }
        }
    }
    
    private static void calculatePotentialParticles() {
        for (int i = 0; i < N; i++) {
            potentialParticles[i][i] = 0;
            for (int j = i + 1; j < N; j++) {
                Particle p1 = particleList.get(i);
                Particle p2 = particleList.get(j);
                double rij = Point.distanceTwo(p1.getR(), p2.getR());
                double potential = eps * (Math.pow(R / rij, 12) - 2 * Math.pow(R / rij, 6));
                potentialParticles[i][j] = potentialParticles[j][i] = potential;
            }
        }
    }
    
    private static void calculatePotentialSum() {
    	calculatePotentialSp();
    	calculatePotentialParticles();
    	potentialSum = 0;
    	for (int i = 1; i < N; i++)
    		for (int j = 0; j < i; j++)
    			potentialSum += potentialParticles[i][j];
    	for (int i = 0; i < N; i++)
    		potentialSum += potentialSp[i];
    	System.out.println(String.format("Suma potencjałów wynosi: %f", potentialSum));
    }

    
    private static void calculatePreasure() {
    	P=0;
    	for(int i=0; i<N; i++) {
    	P+=1/(4*3.14*L*L)*Point.distance(forceSp[i]);
    	}
    	System.out.println(String.format("Ciśnienie na ścianki wynosi: %f Atm", P));
    }
    
    private static void calculateTemperature() {
    	double Esum = 0;
        for (Particle p : particleList) {
            Esum += p.getEnumber();
        }
    	T=2*Esum/(3*N*Particle.kb);
    	System.out.println(String.format("Temperatura wynosi: %f K",T));
    }

    private static void calculateForceParticles() {
        for (int i = 0; i < N; i++) {
            forceParticles[i] = new Point();
            for (int j = 0; j < N; j++) {

                if (j == i)
                    continue;
                else
                {
                Particle p1 = particleList.get(i);
                Particle p2 = particleList.get(j);
                double rij = Point.distanceTwo(p1.getR(), p2.getR());
                double cc = 12 * eps * (Math.pow(R / rij, 12) - Math.pow(R / rij, 6));
                //System.out.println(String.format("cc wynosi %f",cc));
                cc = cc /( rij * rij);
                Point dist = Point.minus(p1.getR(), p2.getR());
                Point force = dist.multiply(cc);
                forceParticles[i] = forceParticles[i].add(force);
                }
            }
        }
    }
    
    private static void calculateForceSp() {
        for (int i = 0; i < N; i++) {
            forceSp[i] = new Point();
            Particle p = particleList.get(i);
            double ri = Point.distance(p.getR());
            if (ri < L) {
            	continue;
            } else {
                double c = f*(L-ri)/ri;
                forceSp[i] = p.getR().multiply(c);
            }
    }
    }
    

    private static void calculateForceSum() {
    	calculateForceParticles();
    	calculateForceSp();
        for(int i=0; i<N; i++) {
            forceSum[i] = forceParticles[i].add(forceSp[i]);
        }
    }

    private static void calculateTotalEnergy() {
    	 Ek = 0;
         for (Particle p : particleList) {
             Ek += p.getEnumber();
         }
         H=Ek+potentialSum;
         System.out.println(String.format("Całkowita Energia wynosi: %f kJ/mol\nHamiltonian wynosi: %f kJ/mol", Ek,H));
    }
    
    private static void checkAverageEnergy() {
        double Esum = 0;
        double Eav=0;
        for (Particle p : particleList) {
            Esum += p.getEnumber();
        }
        Eav = Esum / particleList.size();
        double EavTeor = 3 / 2. * Particle.kb * T;
        System.out.println(String.format("Energia kinetyczna srednia: %f\nTeoria: %f", Eav, EavTeor));
    }

    private static void normalizeMomentum() {
        double pxsum = 0;
        double pysum = 0;
        double pzsum = 0;
        for (Particle p : particleList) {
            pxsum += p.getP().x;
            pysum += p.getP().y;
            pzsum += p.getP().z;
        }

        double pxav = pxsum / particleList.size();
        double pyav = pysum / particleList.size();
        double pzav = pzsum / particleList.size();

        for (Particle p : particleList) {
            Point newP = new Point(p.getP().x - pxav, p.getP().y - pyav, p.getP().z - pzav);
            p.setP(newP);
        }
        System.out.println(String.format("Suma momentów pędu po normalizacji: px=%f, py=%f, pz=%f", pxsum, pysum, pzsum));
    }

    private static void checkMomentum() {
        double pxsum = 0;
        double pysum = 0;
        double pzsum = 0;
        for (Particle p : particleList) {
            pxsum += p.getP().x;
            pysum += p.getP().y;
            pzsum += p.getP().z;
        }
        System.out.println(String.format("Suma momentów pędu: px=%f, py=%f, pz=%f", pxsum, pysum, pzsum));
    }


private static void calculateR() {
    for (Particle p : particleList) {
        Point newPP = new Point(p.getR().x + p.getP().x*tau/m, p.getR().y + p.getP().y*tau/m, p.getR().z + p.getP().z*tau/m);
        p.setR(newPP);     
    	}
    }
	


private static void calculateP() {
	int i=0;
    for (Particle p : particleList) {
    	Point newP = new Point(p.getP().x + forceSum[i].x*0.5*tau, p.getP().y + forceSum[i].y*0.5*tau, p.getP().z + forceSum[i].z*0.5*tau);
        p.setP(newP);
        i++;     
    	}
	}
}

