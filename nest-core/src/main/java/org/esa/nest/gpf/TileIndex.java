package org.esa.nest.gpf;

import org.esa.beam.framework.gpf.Tile;

/**
 * calculates the index into a tile
 */
public final class TileIndex {

    private final int tileOffset;
    private final int tileStride;
    private final int tileMinX;
    private final int tileMinY;

    private int stride = 0;

    public TileIndex(final Tile tile) {
        tileOffset = tile.getScanlineOffset();
        tileStride = tile.getScanlineStride();
        tileMinX = tile.getMinX();
        tileMinY = tile.getMinY();
    }

    public void calculateStride(final int ty) {
        stride = ((ty - tileMinY) * tileStride) + tileOffset;
    }

    public int getIndex(final int tx) {
        return (tx - tileMinX) + stride;
    }
}
