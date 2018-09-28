package edu.mcw.rgd.process;

import java.util.Map;
import java.util.TreeMap;

/**
 * Created by mtutaj on 11/21/2017.
 */
public class Counters {

    private Map<String,Integer> map = new TreeMap<>();

    public void increment(String counter) {
        increment(counter, 1);
    }

    public void increment(String counter, int delta) {
        if( delta!=0 ) {
            int cnt = get(counter);
            map.put(counter, cnt + delta);
        }
    }

    public int get(String counter) {
        Integer cnt = map.get(counter);
        return cnt==null ? 0 : cnt;
    }

    public void dump() {
        for( Map.Entry<String,Integer> entry: map.entrySet() ) {
            System.out.println(entry.getKey()+" : "+entry.getValue());
        }
    }
}
