/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package vrla.util;

import java.util.HashMap;
import java.util.List;

/**
 *
 * @author Sunnyveerla
 */
public class ResultTerm {

    // private float[][] iCorrMatrix;
    //private int[] iRandomVariables;
    private HashMap<String, List> spiralClusters;
    private HashMap<String, List> coreClusters;

    public ResultTerm() {
    }

    /* public void setICorrMatrix(float[][] iCorrMatrix) {
     this.iCorrMatrix = iCorrMatrix;
     }

     public float[][] getICorrMarix() {
     return this.iCorrMatrix;
     }*/

    /*public void setIRandomVariables(int[] iRandomVariables) {
     this.iRandomVariables = iRandomVariables;
     }

     public int[] getIRandomVariables() {
     return this.iRandomVariables;
     }*/
    public void setSpiralClusters(HashMap<String, List> spiralClusters) {
        this.spiralClusters = spiralClusters;
    }

    public HashMap<String, List> getSpiralClusters() {
        return this.spiralClusters;
    }

    public void setCoreClusters(HashMap<String, List> coreClusters) {
        this.coreClusters = coreClusters;
    }

    public HashMap<String, List> getCoreClusters() {
        return this.coreClusters;
    }

    public HashMap<String, List> getClusters(boolean spiral) {
       // System.out.println(spiral);
        if (spiral) {            
            return getSpiralClusters();
        } else {
            return getCoreClusters();
        }
    }
}
