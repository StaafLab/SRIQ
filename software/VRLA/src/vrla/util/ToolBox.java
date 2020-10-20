/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package vrla.util;

import java.awt.Color;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author klin-sve
 */
public class ToolBox {

    public Color[] getDifferentColors() {
        Color[] colorlist = new Color[20];
        colorlist[0] = Color.RED;
        colorlist[1] = Color.BLUE;
        colorlist[2] = Color.GREEN;
        colorlist[3] = Color.ORANGE;
        colorlist[4] = Color.PINK;
        colorlist[5] = Color.MAGENTA;
        colorlist[6] = Color.CYAN;
        colorlist[7] = Color.GRAY;
        colorlist[8] = Color.BLACK;
        colorlist[9] = Color.DARK_GRAY;

        for (int i = 10; i < 20; i++) {
            colorlist[i] = Color.getHSBColor((float) 255 / (i + 1), (float) 255 / ((i + 1) * 2), (float) 255 / ((i + 1) * 3));
        }
        return colorlist;
    }

    public Color[] getDifferentColors(int x) {
        Color[] colorlist = new Color[x];
        colorlist[0] = Color.RED;
        colorlist[1] = Color.BLUE;
        colorlist[2] = Color.GREEN;
        colorlist[3] = Color.ORANGE;
        colorlist[4] = Color.PINK;
        colorlist[5] = Color.MAGENTA;
        colorlist[6] = Color.CYAN;
        colorlist[7] = Color.GRAY;
        colorlist[8] = Color.BLACK;
        colorlist[9] = Color.DARK_GRAY;
        for (int i = 10; i < x; i++) {
            colorlist[i] = Color.getHSBColor((float) 255 / (i + 1), (float) 255 / ((i + 1) * 2), (float) 255 / ((i + 1) * 3));
        }
        return colorlist;
    }

    public String[] getDataFromFile(String inFile) {
        List<String> data = new ArrayList<>();
        try {
            data = Files.readAllLines(Paths.get(inFile));
        } catch (IOException e) {
            e.printStackTrace();
        }
        return data.toArray(new String[0]);
    }

    public double roundDouble(double d, int places) {
        return Math.round(d * Math.pow(10, (double) places)) / Math.pow(10, (double) places);
    }

    public float[][] substract(float results[][], float value) {
//        SStatUtils sutil = new SStatUtils();
        int len = results.length;
        for (int i = 0; i < len; i++) {
            for (int j = i + 1; j < len; j++) {
                results[i][j] = (value - results[i][j]);
                results[j][i] = results[i][j];
            }
        }
        return results;
    }

    public float[][] normalize(float results[][]) {
        int len = results.length;
        float maxD = (float) 0.0;
        for (int i = 0; i < len; i++) {
            for (int j = i + 1; j < len; j++) {
                if (results[i][j] > maxD) {
                    maxD = results[i][j];
                }
            }
        }
        for (int i = 0; i < len; i++) {
            for (int j = i + 1; j < len; j++) {
                //System.out.println(results[i][j]+" DIV "+maxD);
                results[i][j] = (results[i][j] / maxD);
                results[j][i] = results[i][j];
            }
        }
        System.out.println(maxD);
        return results;
    }

    public void writeStringArrayToFile(String data[], String filename) {
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(filename));
            for (int i = 0; i < data.length; i++) {
                if (i == 0) {
                    bw.write(data[i]);
                } else {
                    bw.write("\n" + data[i]);
                }
                bw.flush();
            }
            bw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public String[] getDataFromStringArray(String[] column, int[] row) {
        String[] cell = new String[row.length];
        int count = 0;
        for (int i = 0; i < column.length; i++) {
            if (count == row.length) {
                break;
            }
            for (int j = 0; j < row.length; j++) {
                if (i == row[j]) {
                    cell[j] = column[i];
                    //System.out.println(line+"**********util"+count);
                    count++;
                }
            }

        }
        //System.out.println("Getting cell data");
        return cell;
    }
}
