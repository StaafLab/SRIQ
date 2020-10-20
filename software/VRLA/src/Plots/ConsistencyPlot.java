/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Plots;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import javax.swing.JPanel;
import vrla.util.ToolBox;

/**
 *
 * @author Sunnyveerla
 */
/**
 *
 * @author Sunnyveerla
 */
public class ConsistencyPlot extends JPanel {

    Graphics2D g2d;
    int max;
    Map<String, List<int[]>> cs;

    String title;

    public ConsistencyPlot(Map<String, List<int[]>> cs, int max, String title) {
        this.max = max;
        this.cs = cs;
        this.title = title;
        setSize(2000, 1200);
    }

    public void setTitle(String title) {

        //System.out.println(title);
    }
    public void paint(Graphics g) {
        //System.out.println("Sunny");
        g2d = (Graphics2D) g;
        int x = 20; //distance between the columns
        int tx = 750; //X length
        int ty = 600; //Y length
        int by = 500;
        int size = by / this.max;
        g2d.setColor(Color.BLACK);
        //g2d.drawString("Sunny", x, 50);
        g2d.translate(x, ty - 50);
        g2d.setFont(new Font("SUNNY", Font.BOLD, 12));
        g2d.drawString(this.title, (tx / 2) - this.title.length(), -(by + 10));
        g2d.drawLine(0, 0, 0, -(by + 2 * this.max));
        // g2d.drawLine(0, 0, tx, 0);
        for (int i = 0; i <= this.max; i += 10) {
            g2d.drawString("" + i, -x, -(int) (((float) (i) / this.max) * by));
            //System.out.println(i);
        }
        //Percentage of samples
        int sty = by - 80;
        g2d.drawLine(0, by, 0, 80);
        //g2d.drawLine(0, by, tx, by);
        for (int i = 0; i <= 100; i += 10) {
            g2d.drawString("" + i, -x, by - ((int) (((float) (i) / 100) * sty)));
            //System.out.println(i);
        }
        Object obj[] = cs.keySet().toArray();
        float f[] = new float[obj.length];
        for (int j = 0; j < obj.length; j++) {
            f[j] = Float.parseFloat(obj[j].toString());

        }
        Arrays.sort(f);
        int dy = (size + 2);
        int dx = 21;
        Color c[] = new ToolBox().getDifferentColors();
        for (float o : f) {//System.out.println(o);
            List v = cs.get("" + o);
            int cl = v.size();
            g2d.setColor(Color.BLACK);
            g2d.drawString("" + o, x, 20);
            g2d.drawString("" + o, x, by + 20);
            //System.out.println(o);
            for (int i = 0; i < cl; i++) {
                int[] cp = (int[]) v.get(i);
                if (cp[0] != 0) {
                    int h = (int) ((cp[0] / (float) this.max) * by);
                    g2d.setColor(c[i]);
                    //g2d.setColor(Color.BLACK);
                    g2d.drawRect(x, -h, dx, h);
                    int fill = (int) ((cp[1] * (float) h / 100));
                    g2d.fillRect(x, -fill, dx, fill);
                    //Samples percentage
                    g2d.setColor(Color.BLACK);
                    g2d.drawRect(x, 80, dx, sty);
                    int pfill = (int) ((cp[2] * (float) sty / 100));
                    g2d.fillRect(x, by - pfill, dx, pfill);

                    g2d.drawString("" + (i + 1), (x + 5), -(fill / 2));
                    // g2d.drawString("" + cp[1] + "%", x, -h);
                    x += dx + 5;
                } else {
                    //x += 15;
                }
            }
            x += 20;
        }
        g2d.drawLine(0, 0, x, 0);
        g2d.drawLine(0, by, x, by);      
    }
}
