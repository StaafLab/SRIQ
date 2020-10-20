package vrla.distance;

import java.util.Arrays;

/*
 * JaccardInputFormat.java
 *
 * Created on September 24, 2005, 4:30 AM
 */
/**
 *
 * @author SRIN1
 */
public class JaccardInputFormat {

    /**
     * Creates a new instance of JaccardInputFormat
     */
    public int[][] transformToJaccardInputFormat(float trData[][]) {
        
        int trlen = trData.length;
        int tclen = trData[0].length;
        
        float rData[][] = new float[tclen][trlen];
        for (int x = 0; x < trlen; x++) {
            for (int y = 0; y < tclen; y++) {
                rData[y][x] = trData[x][y];
            }
        }
//System.out.println(trlen + "\t" + tclen);
        int rlen = tclen;
        int clen = trlen;
       //System.out.println(rlen + "\t" + clen);
        int maxNum[] = new int[rlen];
        int totalRows = 0;

        for (int i = 0; i < rlen; i++) {
            int tmp[] = new int[clen];
            for (int j = 0; j < clen; j++) {
                tmp[j] = (int) rData[i][j];
                // System.out.println(tmp[j]);
            }
            Arrays.sort(tmp);
            maxNum[i] = tmp[clen - 1];
            //System.out.println(maxNum[i]);
            totalRows += maxNum[i];
        }
        //System.out.println(totalRows);
        int jData[][] = new int[totalRows][clen];
        int count = 0;
        for (int m = 0; m < rlen; m++) {
            for (int j = 0; j < clen; j++) {
                int tp = (int) rData[m][j];
                if (tp > 0) {
                    //System.out.println(tp);
                    for (int x = 0; x < tp; x++) {
                        jData[count + x][j] = 1;
                        //System.out.print(jData[count+x][j]+"\t");
                    }
                    //System.out.println(count);
                }

            }

            count += maxNum[m];
        }
        /*for (int i = 0; i < totalRows; i++) {
         for (int j = 0; j < clen; j++) {
         System.out.print(jData[i][j] + "\t");
         }
         System.out.println();
         }*/
        return jData;
    }
}
