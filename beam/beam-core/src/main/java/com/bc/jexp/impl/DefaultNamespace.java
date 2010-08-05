/*
 * Copyright (C) 2010 Brockmann Consult GmbH (info@brockmann-consult.de)
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 3 of the License, or (at your option)
 * any later version.
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, see http://www.gnu.org/licenses/
 */

package com.bc.jexp.impl;

import com.bc.jexp.EvalEnv;
import com.bc.jexp.EvalException;
import com.bc.jexp.Namespace;
import com.bc.jexp.Term;

/**
 * A default namespace which includes the constants PI, E and NaN as well as most of the functions from
 * Java {@link Math}, e.g.
 * sin(x),
 * cos(x),
 * tan(x),
 * asin(x),
 * acos(x),
 * atan(x),
 * atan2(x,y),
 * log(x),
 * exp(x),
 * pow(x,y),
 * sqrt(x,y),
 * abs(x),
 * min(x,y),
 * max(x,y),
 * floor(x),
 * rad(x) ({@link Math#toRadians(double)}),
 * deg(x,y) ({@link Math#toDegrees(double)}).
 * inf(x,y) ({@link Double#isInfinite(double)}),
 * nan(x,y) ({@link Double#isNaN(double)}).
 * <p>The following functions are not taken from Java Java {@link Math}:
 * log10(x),
 * exp10(x),
 * sign(x).
 * <p>It also provides the fuzzy comparision functions
 * feq(x, y),
 * fneq(x).
 *
 * @author Norman Fomferra (norman.fomferra@brockmann-consult.de)
 * @version $Revision: 1.3 $ $Date: 2010-08-05 17:00:49 $
 */
public final class DefaultNamespace extends NamespaceImpl {

    private static final double EPS = 1e-6;

    public DefaultNamespace() {
        this(null);
    }

    public DefaultNamespace(Namespace parent) {
        super(parent);
        registerDefaultSymbols();
        registerDefaultFunctions();
    }

    private void registerDefaultSymbols() {
        registerSymbol(SymbolFactory.createConstant("PI", Math.PI));
        registerSymbol(SymbolFactory.createConstant("E", Math.E));
        registerSymbol(SymbolFactory.createConstant("NaN", Double.NaN));
    }

    private void registerDefaultFunctions() {
        registerFunction(new AbstractFunction.D("sin", 1) {

            public double evalD(final EvalEnv env, final Term[] args) {
                return Math.sin(args[0].evalD(env));
            }
        });
        registerFunction(new AbstractFunction.D("cos", 1) {

            public double evalD(final EvalEnv env, final Term[] args) {
                return Math.cos(args[0].evalD(env));
            }
        });
        registerFunction(new AbstractFunction.D("tan", 1) {

            public double evalD(final EvalEnv env, final Term[] args) {
                return Math.tan(args[0].evalD(env));
            }
        });
        registerFunction(new AbstractFunction.D("asin", 1) {

            public double evalD(final EvalEnv env, final Term[] args) {
                return Math.asin(args[0].evalD(env));
            }
        });
        registerFunction(new AbstractFunction.D("acos", 1) {

            public double evalD(final EvalEnv env, final Term[] args) {
                return Math.acos(args[0].evalD(env));
            }
        });
        registerFunction(new AbstractFunction.D("atan", 1) {

            public double evalD(final EvalEnv env, final Term[] args) {
                return Math.atan(args[0].evalD(env));
            }
        });
        registerFunction(new AbstractFunction.D("atan2", 2) {

            public double evalD(final EvalEnv env, final Term[] args) {
                return Math.atan2(args[0].evalD(env), args[1].evalD(env));
            }
        });
        registerFunction(new AbstractFunction.D("log", 1) {

            public double evalD(final EvalEnv env, final Term[] args) {
                return Math.log(args[0].evalD(env));
            }
        });

        registerFunction(new AbstractFunction.D("log10", 1) {
            public double evalD(EvalEnv env, Term[] args) throws EvalException {
                return Math.log10(args[0].evalD(env));
            }
        });

        registerFunction(new AbstractFunction.D("exp", 1) {

            public double evalD(final EvalEnv env, final Term[] args) {
                return Math.exp(args[0].evalD(env));
            }
        });

        registerFunction(new AbstractFunction.D("exp10", 1) {
            public double evalD(EvalEnv env, Term[] args) throws EvalException {
                return Math.pow(10.0, args[0].evalD(env));
            }
        });

        registerFunction(new AbstractFunction.D("sqrt", 1) {

            public double evalD(final EvalEnv env, final Term[] args) {
                return Math.sqrt(args[0].evalD(env));
            }
        });
        registerFunction(new AbstractFunction.D("pow", 2) {

            public double evalD(final EvalEnv env, final Term[] args) {
                return Math.pow(args[0].evalD(env), args[1].evalD(env));
            }
        });
        registerFunction(new AbstractFunction.I("min", 2) {

            public int evalI(final EvalEnv env, final Term[] args) {
                return Math.min(args[0].evalI(env), args[1].evalI(env));
            }
        });
        registerFunction(new AbstractFunction.D("min", 2) {

            public double evalD(final EvalEnv env, final Term[] args) {
                return Math.min(args[0].evalD(env), args[1].evalD(env));
            }
        });
        registerFunction(new AbstractFunction.I("max", 2) {

            public int evalI(final EvalEnv env, final Term[] args) {
                return Math.max(args[0].evalI(env), args[1].evalI(env));
            }
        });
        registerFunction(new AbstractFunction.D("max", 2) {

            public double evalD(final EvalEnv env, final Term[] args) {
                return Math.max(args[0].evalD(env), args[1].evalD(env));
            }
        });
        registerFunction(new AbstractFunction.D("floor", 1) {
            public double evalD(EvalEnv env, Term[] args) throws EvalException {
                return Math.floor(args[0].evalD(env));
            }
        });
        registerFunction(new AbstractFunction.D("round", 1) {
            public double evalD(EvalEnv env, Term[] args) throws EvalException {
                return Math.round(args[0].evalD(env));
            }
        });
        registerFunction(new AbstractFunction.D("ceil", 1) {
            public double evalD(EvalEnv env, Term[] args) throws EvalException {
                return Math.ceil(args[0].evalD(env));
            }
        });
        registerFunction(new AbstractFunction.D("rint", 1) {
            public double evalD(EvalEnv env, Term[] args) throws EvalException {
                return Math.rint(args[0].evalD(env));
            }
        });
        registerFunction(new AbstractFunction.I("sign", 1) {

            public int evalI(final EvalEnv env, final Term[] args) {
                return ExtMath.sign(args[0].evalI(env));
            }
        });
        registerFunction(new AbstractFunction.D("sign", 1) {

            public double evalD(final EvalEnv env, final Term[] args) {
                return ExtMath.sign(args[0].evalD(env));
            }
        });
        registerFunction(new AbstractFunction.I("abs", 1) {

            public int evalI(final EvalEnv env, final Term[] args) {
                return Math.abs(args[0].evalI(env));
            }
        });
        registerFunction(new AbstractFunction.D("abs", 1) {

            public double evalD(final EvalEnv env, final Term[] args) {
                return Math.abs(args[0].evalD(env));
            }
        });
        registerFunction(new AbstractFunction.D("deg", 1) {

            public double evalD(final EvalEnv env, final Term[] args) {
                return Math.toDegrees(args[0].evalD(env));
            }
        });
        registerFunction(new AbstractFunction.D("rad", 1) {

            public double evalD(final EvalEnv env, final Term[] args) {
                return Math.toRadians(args[0].evalD(env));
            }
        });
        registerFunction(new AbstractFunction.D("ampl", 2) {

            public double evalD(final EvalEnv env, final Term[] args) {
                final double a = args[0].evalD(env);
                final double b = args[1].evalD(env);
                return Math.sqrt(a * a + b * b);
            }
        });
        registerFunction(new AbstractFunction.D("phase", 2) {

            public double evalD(final EvalEnv env, final Term[] args) {
                final double a = args[0].evalD(env);
                final double b = args[1].evalD(env);
                return Math.atan2(b, a);
            }
        });
        registerFunction(new AbstractFunction.B("feq", 2) {
            public boolean evalB(EvalEnv env, Term[] args) throws EvalException {
                final double x1 = args[0].evalD(env);
                final double x2 = args[1].evalD(env);
                return ExtMath.feq(x1, x2, EPS);
            }
        });
        registerFunction(new AbstractFunction.B("feq", 3) {
            public boolean evalB(EvalEnv env, Term[] args) throws EvalException {
                final double x1 = args[0].evalD(env);
                final double x2 = args[1].evalD(env);
                final double eps = args[2].evalD(env);
                return ExtMath.feq(x1, x2, eps);
            }
        });
        registerFunction(new AbstractFunction.B("fneq", 2) {
            public boolean evalB(EvalEnv env, Term[] args) throws EvalException {
                final double x1 = args[0].evalD(env);
                final double x2 = args[1].evalD(env);
                return ExtMath.fneq(x1, x2, EPS);
            }
        });
        registerFunction(new AbstractFunction.B("fneq", 3) {
            public boolean evalB(EvalEnv env, Term[] args) throws EvalException {
                final double x1 = args[0].evalD(env);
                final double x2 = args[1].evalD(env);
                final double eps = args[2].evalD(env);
                return ExtMath.fneq(x1, x2, eps);
            }
        });
        registerFunction(new AbstractFunction.B("inf", 1) {
            public boolean evalB(EvalEnv env, Term[] args) throws EvalException {
                return Double.isInfinite(args[0].evalD(env));
            }
        });
        registerFunction(new AbstractFunction.B("nan", 1) {
            public boolean evalB(EvalEnv env, Term[] args) throws EvalException {
                return Double.isNaN(args[0].evalD(env));
            }
        });
    }

}
