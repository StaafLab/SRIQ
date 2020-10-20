package vrla.distance;

/*
 * STM.java
 *
 * Created on April 12, 2005, 5:24 PM
 */
/**
 *
 * @author SRIN1
 */
// This program is to create the Jacknife Correlation Matrix. Using the output file from this program as input to QTClust program to get QTC clusters.

public class JacknifeCorrelationDistanceMatrix {

    public float[][] getJacknifeCorrelationDistanceMatrix(float fdata[][]) {
        int rlen = fdata.length;
        int clen = fdata[0].length;        
        float cdata[][] = new float[rlen][rlen];
        float r = (float) 0.0;
        for (int i = 0; i < rlen; i++) {
            cdata[i][i] = 0;
            for (int o = i + 1; o < rlen; o++) {              
                float sAvg = (float) 0.0;
                for (int k = -1; k < clen; k++) {
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
                    sAvg += r;                   

                }                
                cdata[i][o] = (1 - (sAvg / clen));//convert to distance 1-r
                cdata[o][i] = cdata[i][o];
            }
        }
        return cdata;
    } 
}
