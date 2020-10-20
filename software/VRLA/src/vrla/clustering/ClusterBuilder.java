/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package vrla.clustering;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import vrla.input.PropertiesConfiguration;

/**
 *
 * @author klin-sve
 */
public class ClusterBuilder {
    public synchronized static HashMap<String, HashMap<String, List<List<Integer>>>> getQTClusters(float sData[][], double cutOffs[]) {
        HashMap<String, HashMap<String, List<List<Integer>>>> ht = new HashMap<>();
        HashMap<String, List<List<Integer>>> ch = new HashMap<>(cutOffs.length);
        HashMap<String, List<List<Integer>>> sh = new HashMap<>(cutOffs.length);
        int num_of_items = sData.length;
        for (double d : cutOffs) {
           // STC stc = new STC(sData, d, cSize);
            //SQTClustering stc = new SQTClustering(sData, d, cSize);
//            List<Integer> allUniqueIDIndices = new ArrayList<>();
//            allUniqueIDIndices = Collections.synchronizedList(allUniqueIDIndices);
            List<List<Integer>> clusterArrayList = new ArrayList<>();

//            for (int i = 0; i < num_of_items; i++) {
//                allUniqueIDIndices.add(i);
//            }
            clusterArrayList = new STClusterAnalysis(sData, d, PropertiesConfiguration.SAMPLE_BAG_SIZE).getClusterAnalysisClassObject().getAllClusters();
            ch.put("" + d, clusterArrayList);
            int cvs = clusterArrayList.size();
            if (PropertiesConfiguration.SPIRAL && !clusterArrayList.get(cvs - 1).isEmpty()) {
                sh.put("" + d, getSpiralClusters(clusterArrayList, sData));
            }
        }
        ht.put("ch", ch);
        ht.put("sh", sh);
        return ht;
    }
    

    public synchronized static List<List<Integer>> getSpiralClusters(List<List<Integer>> clusterArrayList, float[][] fdata) {
        int cSize = clusterArrayList.size();
        List<Integer> unassigned = (ArrayList) clusterArrayList.get(cSize - 1);
        int ulen = unassigned.size();
        if (ulen == 0) {
            return clusterArrayList;
        }
        HashMap<String, List> aCluster = new HashMap(cSize - 1);

        for (int i = 0; i < ulen; i++) {
            float avgDistance = Float.POSITIVE_INFINITY;
            int addToCluster = 0;
            int u = unassigned.get(i);
            for (int j = 0; j < cSize - 1; j++) {
                Object sp[] = clusterArrayList.get(j).toArray();
                int len = sp.length;
                float tmp = (float) 0.0;
                for (int o = 0; o < len; o++) {
                    tmp += fdata[u][Integer.parseInt(sp[o].toString())];
                }
                float avgTmp = (tmp / len);
                if (avgDistance > avgTmp) {
                    avgDistance = avgTmp;
                    //System.out.println(avgDistance);
                    addToCluster = j; //j starts from zero
                }
            }
            if (!aCluster.containsKey("" + addToCluster)) {
                ArrayList vTmp = new ArrayList();
                vTmp.add(u);
                aCluster.put("" + addToCluster, vTmp);
            } else {
                List vTmp = aCluster.get("" + addToCluster);
                vTmp.add(u);
                aCluster.put("" + addToCluster, vTmp);
            }
        }

        System.out.println("Total Clusters " + (cSize - 1));
        List<List<Integer>> spiralClusters = new ArrayList<>(cSize);
        for (int x = 0; x < cSize - 1; x++) {
            List addToArrayList = new ArrayList(clusterArrayList.get(x));
            if (aCluster.containsKey("" + x)) {
                addToArrayList.addAll(aCluster.get("" + x));
            }
            spiralClusters.add(x, addToArrayList);
        }
        spiralClusters.add((cSize - 1), unassigned);
        return (spiralClusters);
    }
}


/*private ResultTerm getQTClusters(float sData[][], float distCutOffs[], int cSize, boolean spiral) {
     ResultTerm rt = new ResultTerm();
     HashMap<String, ArrayList> ch = new HashMap(distCutOffs.length);
     HashMap<String, ArrayList> sh = new HashMap(distCutOffs.length);
     int num_of_items = sData.length;
     for (float d : distCutOffs) {
     STC stc = new STC(sData, d, cSize);
     ArrayList allUniqueIDIndices = new ArrayList();
     ArrayList<ArrayList<Integer>> clusterArrayList = new ArrayList();

     for (int i = 0; i < num_of_items; i++) {
     allUniqueIDIndices.add(new Integer(i));
     }
     clusterArrayList = stc.getAllClusters(allUniqueIDIndices);
     ch.put("" + d, clusterArrayList);
     int cvs = clusterArrayList.size();
     if (spiral && clusterArrayList.get(cvs - 1).size() != 0) {
     sh.put("" + d, getSpiralClusters(clusterArrayList, sData));
     }
     }
     rt.setCoreClusters(ch);
     rt.setSpiralClusters(sh);
     return rt;
     }*/
