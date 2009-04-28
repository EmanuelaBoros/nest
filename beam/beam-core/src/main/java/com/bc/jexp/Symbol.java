/*
 * $Id: Symbol.java,v 1.1 2009-04-28 14:39:32 lveci Exp $
 *
 * Copyright (C) 2002 by Brockmann Consult (info@brockmann-consult.de)
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

package com.bc.jexp;


/**
 * Represents a read-only symbol. A symbol can be a named constant or variable.
 * It has a return type an can be evaluated.
 *
 * <p>Within an expression, a reference to a symbol is created if the parser
 * encounters a name and this name can be resolved through the parser's current namespace.
 * The resulting term in this case is an instance of <code>{@link Term.Ref}</code>.
 * @author Norman Fomferra (norman.fomferra@brockmann-consult.de)
 * @version $Revision: 1.1 $ $Date: 2009-04-28 14:39:32 $
 */
public interface Symbol {

    /**
     * Gets the symbol's name.
     * @return the name, should never be <code>null</code>
     */
    String getName();

    /**
     * Gets the symbol's return type.
     * @return the type, should always be one of the <code>TYPE_</code>X constants
     *         defined in the <code>Term</code> class.
     */
    int getRetType();

    /**
     * Evaluates this symbol to a <code>boolean</code> value.
     * @param env the application dependant environment.
     * @return a <code>boolean</code> value
     * @throws EvalException if the evaluation fails
     */
    boolean evalB(EvalEnv env) throws EvalException;

    /**
     * Evaluates this symbol to an <code>int</code> value.
     * @param env the application dependant environment.
     * @return an <code>int</code> value
     * @throws EvalException if the evaluation fails
     */
    int evalI(EvalEnv env) throws EvalException;

    /**
     * Evaluates this symbol to a <code>double</code> value.
     * @param env the application dependant environment.
     * @return a <code>double</code> value
     * @throws EvalException if the evaluation fails
     */
    double evalD(EvalEnv env) throws EvalException;

    /**
     * Evaluates this symbol to a <code>String</code> value.
     * @param env the application dependant environment.
     * @return a <code>double</code> value
     * @throws EvalException if the evaluation fails
     */
    String evalS(EvalEnv env) throws EvalException;
}
