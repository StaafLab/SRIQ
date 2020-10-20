/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package vrla.distance;

/**
 *
 * @author sunnyveerla
 */
public class SquaredEucledianDistanceMatrix {

    public SquaredEucledianDistanceMatrix() {

    }

    public float[][] getEucledianDistanceMatrix(float fdata[][]) {
        int rlen = fdata.length;
        int clen = fdata[0].length;
        float cdata[][] = new float[rlen][rlen];
        float d = Float.POSITIVE_INFINITY;
        for (int i = 0; i < rlen; i++) {
            cdata[i][i] = 0;
            for (int o = i + 1; o < rlen; o++) {
                //for (int k = -1; k < ncol - 1; k++) {
                float SUMX = (float) 0.0;
                float SUMY = (float) 0.0;
                float SUMXY = (float) 0.0;
                float SUMX2 = (float) 0.0;
                float SUMY2 = (float) 0.0;
                float sdx = (float) 0.0;
                float sdy = (float) 0.0;
                for (int j = 0; j < clen; j++) {
                    SUMX += fdata[i][j];
                    SUMX2 += fdata[i][j] * fdata[i][j];
                    SUMY += fdata[o][j];
                    SUMY2 += fdata[o][j] * fdata[o][j];
                    SUMXY += fdata[i][j] * fdata[o][j];
                }
               
                 d = (float) Math.sqrt(SUMX2 + SUMY2 - (2 * SUMXY));
                // System.out.println(d);
                    
                    if (("" + d).equals("NaN")) {
                        d = Float.POSITIVE_INFINITY;
                    }
                cdata[i][o] = d;
                cdata[o][i] = d;
            }
        }
        return cdata;
    }

    /* public static void main(String[] args) {
     String inpath = "/Users/sunnyveerla/Sunny/Projects/Illumina/OncoChip/Data Analysis/";
     String filename = "Promoter_SNP_data";
     String infile = inpath + filename + ".txt";
     //String tfile = "/users/sunnyveerla/sunnytf.txt";

     UtilInterface util = new UtilInterface();
     int nrow = util.lineCount(infile);
     int ncol = util.colCount();
     int count = 0;
     BufferedReader br = null;
     BufferedWriter bw = null;
     BufferedWriter bwc = null;
     String line = null;
     String genes[] = new String[nrow - 1];
     float fdata[][] = new float[nrow - 1][ncol - 1];
     int lc = 0;

     float nc = (float) ncol - 1;
     float r = (float) 0.0;
     float tr = (float) 0.5;
     String tab = "";
     String outfile = inpath + filename + "_QTC.txt";
     String coutfile = inpath + filename + "_QTC_cytoscape.txt";
     try {
     br = new BufferedReader(new FileReader(infile));
     bw = new BufferedWriter(new FileWriter(outfile));
     bwc = new BufferedWriter(new FileWriter(coutfile));
     while ((line = br.readLine()) != null) {
     if (lc == 0) {
     lc++;
     continue;
     }
     System.out.println(count);
     String sp[] = line.split("\t");
     genes[count] = sp[0];
     for (int i = 1; i < ncol; i++) {
     fdata[count][i - 1] = Float.parseFloat(sp[i]);
     }
     count++;
     sp = null;
     }
     for (int i = 0; i < nrow - 1; i++) {
     bw.write("1");
     bw.flush();
     for (int o = i + 1; o < nrow - 1; o++) {

     //for (int k = -1; k < ncol - 1; k++) {
     float SUMX = (float) 0.0;
     float SUMY = (float) 0.0;
     float SUMXY = (float) 0.0;
     float SUMX2 = (float) 0.0;
     float SUMY2 = (float) 0.0;
     float sdx = (float) 0.0;
     float sdy = (float) 0.0;
     for (int j = 0; j < ncol - 1; j++) {
     /*if (k != -1 && j == k) {
     continue;
     }*/
    /*SUMX += fdata[i][j];
     SUMX2 += fdata[i][j] * fdata[i][j];
     SUMY += fdata[o][j];
     SUMY2 += fdata[o][j] * fdata[o][j];
     SUMXY += fdata[i][j] * fdata[o][j];
     }
     sdx = (nc * SUMX2 - (float) Math.pow(SUMX, 2));
     sdy = (nc * SUMY2 - (float) Math.pow(SUMY, 2));
     r = (nc * SUMXY - SUMX * SUMY) / (float) Math.sqrt(sdx * sdy);
     //System.out.println(r + "XXXXXXX" + i + "*******" + o);
     if (("" + r).equals("NaN") && sdx == sdy) {
     r = (float) 1.0;
     } else if (("" + r).equals("NaN") && sdx != sdy) {
     r = (float) 0.0;
     }

     //}
     //System.out.print(tr + "\t");
     if (r >= tr) {
     bwc.write(genes[i] + "\t" + genes[o] + "\t" + (1 - r) + "\n");
     bwc.flush();
     }
     bw.write("\t" + r);
     bw.flush();
     }
     System.out.println(i + "****");
     if (i < nrow - 2) {
     tab += "\t";
     bw.write("\n" + tab);
     }
     bw.flush();
     }
     bwc.close();
     bw.close();
     } catch (IOException e) {
     e.printStackTrace();
     }

     }*/

    
}
