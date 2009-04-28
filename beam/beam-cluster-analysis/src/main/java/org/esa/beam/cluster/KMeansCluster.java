/* 
 * Copyright (C) 2002-2008 by Brockmann Consult
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation. This program is distributed in the hope it will
 * be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */
package org.esa.beam.cluster;

/**
 * todo - add API doc
 *
 * @author Ralf Quast
 * @version $Revision: 1.1 $ $Date: 2009-04-28 14:37:13 $
 * @since BEAM 4.2
 */
class KMeansCluster {
    private final double[] mean;
    private final int memberCount;

    KMeansCluster(double[] mean, int memberCount) {
        this.mean = mean;
        this.memberCount = memberCount;
    }

    double[] getMean() {
        return mean;
    }

    int getMemberCount() {
        return memberCount;
    }
}
