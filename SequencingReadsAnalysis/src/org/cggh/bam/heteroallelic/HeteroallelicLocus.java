package org.cggh.bam.heteroallelic;

import org.cggh.common.genome.*;

public class HeteroallelicLocus {
	
    private String       name;
    private GenomeRegion region;
    private int          startCodon;
    private boolean      isReverse;

    public HeteroallelicLocus(String name, GenomeRegion region, int startCodon, boolean isReverse) {
        this.name = name;
        this.region = region;
        this.startCodon = startCodon;
        this.isReverse = isReverse;
    }

    public String getName() {
        return this.name;
    }

    public GenomeRegion getRegion() {
        return this.region;
    }

    public int getStartCodon() {
        return this.startCodon;
    }

    public boolean isReverse() {
        return this.isReverse;
    }
}

