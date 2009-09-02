package com.bc.ceres.binio.internal;

import com.bc.ceres.binio.*;

import java.io.IOException;


final class SimpleMember implements MemberInstance {
    private final DataContext context;
    private final CollectionData parent;
    private final SimpleType simpleType;
    private final Segment segment;
    private final int segmentOffset;
    private DataAccessor dataAccessor;

    protected SimpleMember(DataContext context,
                           CollectionData parent,
                           SimpleType simpleType,
                           Segment segment,
                           int segmentOffset) {
        this.context = context;
        this.parent = parent;
        this.simpleType = simpleType;
        this.segment = segment;
        this.segmentOffset = segmentOffset;
    }

    public CollectionData getParent() {
        return parent;
    }

    public Segment getSegment() {
        return segment;
    }

    public int getSegmentOffset() {
        return segmentOffset;
    }

    @Override
    public long getPosition() {
        return segment.getPosition() + segmentOffset;
    }

    @Override
    public long getSize() {
        return simpleType.getSize();
    }

    @Override
    public boolean isSizeResolved() {
        return true;
    }

    @Override
    public void resolveSize() {
        // ok
    }

    @Override
    public Type getType() {
        return simpleType;
    }

    ////////////////////////////////////////////////////
    // data access

    @Override
    public byte getByte() throws IOException {
        ensureDataAccessible();
        return dataAccessor.getByte(segment.getData(), segmentOffset);
    }

    @Override
    public void setByte(byte value) throws IOException {
        ensureDataAccessible();
        dataAccessor.setByte(segment.getData(), segmentOffset, value);
        segment.setDirty(true);
    }

    @Override
    public short getShort() throws IOException {
        ensureDataAccessible();
        return dataAccessor.getShort(segment.getData(), segmentOffset);
    }

    @Override
    public void setShort(short value) throws IOException {
        ensureDataAccessible();
        dataAccessor.setShort(segment.getData(), segmentOffset, value);
        segment.setDirty(true);
    }

    @Override
    public int getInt() throws IOException {
        ensureDataAccessible();
        return dataAccessor.getInt(segment.getData(), segmentOffset);
    }

    @Override
    public void setInt(int value) throws IOException {
        ensureDataAccessible();
        dataAccessor.setInt(segment.getData(), segmentOffset, value);
        segment.setDirty(true);
    }

    @Override
    public long getLong() throws IOException {
        ensureDataAccessible();
        return dataAccessor.getLong(segment.getData(), segmentOffset);
    }

    @Override
    public void setLong(long value) throws IOException {
        ensureDataAccessible();
        dataAccessor.setLong(segment.getData(), segmentOffset, value);
        segment.setDirty(true);
    }

    @Override
    public float getFloat() throws IOException {
        ensureDataAccessible();
        return dataAccessor.getFloat(segment.getData(), segmentOffset);
    }

    @Override
    public void setFloat(float value) throws IOException {
        ensureDataAccessible();
        dataAccessor.setFloat(segment.getData(), segmentOffset, value);
        segment.setDirty(true);
    }

    @Override
    public double getDouble() throws IOException {
        ensureDataAccessible();
        return dataAccessor.getDouble(segment.getData(), segmentOffset);
    }

    @Override
    public void setDouble(double value) throws IOException {
        ensureDataAccessible();
        dataAccessor.setDouble(segment.getData(), segmentOffset, value);
        segment.setDirty(true);
    }

    @Override
    public SequenceInstance getSequence() {
        throw new DataAccessException();
    }

    @Override
    public CompoundInstance getCompound() {
        throw new DataAccessException();
    }

    @Override
    public void flush() throws IOException {
        segment.flushData(context);
    }

    // data access
    ////////////////////////////////////////////////////

    private void ensureDataAccessible() throws IOException {
        if (dataAccessor == null) {
            this.dataAccessor = DataAccessor.getInstance(simpleType, context.getFormat().getByteOrder());
            if (!segment.isDataAccessible()) {
                segment.makeDataAccessible(context);
            }
        }
    }
}
