/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package vrla.input;

import java.io.FileReader;
import java.io.IOException;

/**
 *
 * @author klin-sve
 */
public class LoadProperties {

    public static LinkedProperties PROP;
    public String propFilePath;

    public LoadProperties(String propFilePath) {
        this.propFilePath = propFilePath;
        setPropertiesConfiguration();
    }

    private void setPropertiesConfiguration() {
        PROP = new LinkedProperties();
        try {
            PROP.load(new FileReader(propFilePath));
        } catch (IOException e) {
            System.out.println("IOException in LoadProperties");
        }
    }

}
