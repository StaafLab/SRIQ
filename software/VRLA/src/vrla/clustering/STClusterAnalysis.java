/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package vrla.clustering;

import vrla.clustering.SQTC;

/**
 *
 * @author Sunnyveerla
 */
public class STClusterAnalysis {

    private float fdata[][];
    private double adjustedDiameter;
    private int minimumClusterSize;
    private int number_of_items;

    public STClusterAnalysis(float fdata[][], double adjustedDiameter, int minimumClusterSize) {
        this.adjustedDiameter = adjustedDiameter;
        this.minimumClusterSize = minimumClusterSize;
        this.fdata = fdata;
        this.number_of_items = fdata.length;
    }

    public synchronized SQTC getClusterAnalysisClassObject() {
        //return new QTC(fdata, adjustedDiameter, minimumClusterSize);
        return new SQTC(fdata, adjustedDiameter, minimumClusterSize);
    }
}
