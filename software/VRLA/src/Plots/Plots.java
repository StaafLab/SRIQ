package Plots;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.GraphicsConfiguration;
import java.awt.GraphicsDevice;
import java.awt.GraphicsEnvironment;
import java.awt.Transparency;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import javax.imageio.ImageIO;
import javax.swing.JPanel;
import vrla.util.ToolBox;
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author sunnyveerla
 */
public class Plots {

    /**
     * @param args the command line arguments
     */
//    UtilInterface util = new UtilInterface();
    String Header = "";
    JPanel plot;

    public Plots() {
        plot = new JPanel();
        //plot.setSize(1000, 2000);
        plot.setBackground(Color.WHITE);
    }

    public void makeConsistencyPlot(String infileA, String title) {
        String data[] = new ToolBox().getDataFromFile(infileA);
        int len = data.length;
        Map<String, List<int[]>> hc = new HashMap(len - 1);
        String st[] = data[0].split("\t");
        System.out.println(st.length);
        int max = 0;
        for (int i = 1; i < st.length; i++) {
            List cv = new ArrayList(len - 1);
            for (int j = 1; j < len; j++) {
                String sp[] = data[j].split("\t");
                String c[] = sp[i].split("_");

                int v[] = new int[3];
                if (c.length > 1) {
                    //c[1] = c[1].substring(0, c[1].length() - 1);
                    v[0] = Integer.parseInt(c[0]);
                    if (i == 1) {
                        max += v[0];
                    }
                    v[1] = (int) (Float.parseFloat(c[2]) * 100);
                    v[2] = (int) (Float.parseFloat(c[4]) * 100);
                }
                cv.add((j - 1), v);
            }
            hc.put("" + Float.parseFloat(st[i].trim()), cv);
        }
        System.out.println("Its here");
        ConsistencyPlot cp = new ConsistencyPlot(hc, max, title);
        cp.setPreferredSize(new Dimension(1500, 1200));
        plot.setSize(new Dimension(2000, 1300));
        plot.add(cp);
        
        saveComponentAsPNG(plot, infileA.substring(0, infileA.indexOf(".txt")));
        //display();
    }

    public void saveComponentAsPNG(JPanel plots, String filename) {
        //MessageBox msg = new MessageBox();
        //msg.showMessage("Image saving in process", 1);
        GraphicsEnvironment ge = GraphicsEnvironment.getLocalGraphicsEnvironment();
        GraphicsDevice gs = ge.getDefaultScreenDevice();
        GraphicsConfiguration gc = gs.getDefaultConfiguration();
        //BufferedImage myImage = new BufferedImage(getWidth(), getHeight(),BufferedImage.TYPE_4BYTE_ABGR);
        BufferedImage myImage = gc.createCompatibleImage(plots.getWidth(), plots.getHeight(), Transparency.TRANSLUCENT);
        Graphics2D g2 = myImage.createGraphics();

        plots.paint(g2);
        g2.dispose();
        try {
            ImageIO.write(myImage, "png", new File(filename + ".png"));
        } catch (Exception e) {
            System.out.println(e);
        }
        //msg.hide();
        // msg.showMessage("Image saving process is finished", 1);
    }

}
