/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package vrla.util;

import java.util.ArrayList;
import java.util.Random;
import vrla.data.FileParser;
import vrla.input.PropertiesConfiguration;

/**
 *
 * @author klin-sve
 */
public class BaggingProcess {

    private static int[] getRandomNumbers() {
        Random rndm = new Random();
        int[] rn = new int[PropertiesConfiguration.GENE_BAG_SIZE];
        ArrayList<Integer> rv = new ArrayList<>(PropertiesConfiguration.GENE_BAG_SIZE);
        int i = 0;
        while (rv.size() < PropertiesConfiguration.GENE_BAG_SIZE) {
            int tmp = rndm.nextInt(PropertiesConfiguration.TOTAL_SAMPLES);
            if (!rv.contains(tmp)) {
                rn[i] = tmp;
                rv.add(tmp);
                i++;
            }
        }
        return rn;
    }

    public static float[][] getDataBag() {
        float rdata[][] = new float[PropertiesConfiguration.TOTAL_SAMPLES][PropertiesConfiguration.GENE_BAG_SIZE];
        int[] rn = getRandomNumbers();
        for (int g = 0; g < PropertiesConfiguration.GENE_BAG_SIZE; g++) {
            //System.out.println(rm);
            float tmp[] = FileParser.hData.get(rn[g]);
            for (int s = 0; s < PropertiesConfiguration.TOTAL_SAMPLES; s++) {
                rdata[s][g] = tmp[s];
            }
        }
        return rdata;
    }
}

