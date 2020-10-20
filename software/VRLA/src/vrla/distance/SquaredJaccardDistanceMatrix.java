package vrla.distance;

/*
 * SquaredJaccardDistanceMatrix.java
 *
 * Created on September 24, 2005, 4:24 AM
 */
/**
 *
 * @author SRIN1
 */
public class SquaredJaccardDistanceMatrix {    
    public float[][] getJaccardDistance(float fData[][]) {
        int[][] jData = new JaccardInputFormat().transformToJaccardInputFormat(fData);
        int row = jData.length;
        int col = jData[0].length;
        float[][] jDMatrix = new float[col][col];
        int count = 0;
        int a = 0;
        int b = 0;
        int c = 0;
        for (int j = 0; j < col; j++) {

            //System.out.println(hdr[j + 1]);
            for (int k = j; k < col; k++) {
                if (j == k) {
                    jDMatrix[j][k] = 0;
                    continue;
                }
                for (int i = 0; i < row; i++) {
                    if (jData[i][j] == 1 && jData[i][k] == 1) {
                        a += 1;
                    } else if (jData[i][j] == 1 && jData[i][k] == 0) {
                        b += 1;
                    } else if (jData[i][j] == 0 && jData[i][k] == 1) {
                        c += 1;
                    }
                }
                jDMatrix[j][k] = (1 - (a / ((float) a + b + c)));
                jDMatrix[k][j] = jDMatrix[j][k];
                a = 0;
                b = 0;
                c = 0;
            }
        }
        /*for (int i = 0; i < col; i++) {
            for (int j = 0; j < col; j++) {
                System.out.print(jDMatrix[i][j]+"\t");
            }
            System.out.println();
        }*/
        return jDMatrix;
    }

}
