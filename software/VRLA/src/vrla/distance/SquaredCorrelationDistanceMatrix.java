/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package vrla.distance;

/**
 *
 * @author sunnyveerla
 */
public class SquaredCorrelationDistanceMatrix {

    public synchronized float[][] getCorrelationDistanceMatrix(float fdata[][]) {
//        System.out.println("In Correlation distance matrix");
        int rlen = fdata.length;
        int clen = fdata[0].length;
        float cdata[][] = new float[rlen][rlen];
        float r = (float) 0.0;
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
                sdx = (clen * SUMX2 - (float) Math.pow(SUMX, 2));
                sdy = (clen * SUMY2 - (float) Math.pow(SUMY, 2));
                //System.out.println(sdx + "\t" + sdy);
                r = (clen * SUMXY - SUMX * SUMY) / (float) Math.sqrt(sdx * sdy);
                if (("" + r).equals("NaN") && sdx == sdy) {
                    r = (float) 1.0;
                } else if (("" + r).equals("NaN") && sdx != sdy) {

                    r = (float) 0.0;

                }

                //float tmpR = (float)sutil.roundDouble((double)(1-r),3);
                cdata[i][o] = (1 - r);//convert to distance 1-r
                cdata[o][i] = (1 - r);
            }
        }
        return cdata;
    }
}
