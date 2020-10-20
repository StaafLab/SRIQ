/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package vrla.data;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 *
 * @author klin-sve
 */
public class FileParser {

    public static String filePath;
    public static String[] header;
    public static String[] rownames;
    public static HashMap<Integer, float[]> hData;

    /**
     * @return the hData
     */
    public static HashMap<Integer, float[]> gethData() {
        return hData;
    }

    public FileParser(String filePath) {
        this.filePath = filePath;
        this.hData = new HashMap<>();
        sethData(hData);
    }

    /**
     * @param ahData the hData to set
     */
    public void sethData(HashMap<Integer, float[]> ahData) {
        String line = null;
        int h = 0;
        int lc = 0;
        List<String> rows = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(filePath))) {
            while ((line = br.readLine()) != null) {
                if (h == 0) {
//                    System.out.println(line);
//                    line = line.substring(line.indexOf("\t"), line.length()).trim();
                    header = line.split("\t");
                    h++;
                    continue;
                }
                String[] sp = line.split("\t");
                float[] fdata = new float[header.length];
                for (int i = 0; i < sp.length - 1; i++) {
                    fdata[i] = Float.parseFloat(sp[i + 1]);
                }
//                rownames[lc] = sp[0];
                rows.add(sp[0]);
                hData.put(lc++, fdata);
            }
            rownames = rows.stream().toArray(String[]::new);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
