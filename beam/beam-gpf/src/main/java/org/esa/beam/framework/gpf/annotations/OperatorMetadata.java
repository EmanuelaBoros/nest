package org.esa.beam.framework.gpf.annotations;

import java.lang.annotation.*;

/**
 * Provides metadata for an operator. This annotation is valid for
 * extensions of the {@link org.esa.beam.framework.gpf.Operator Operator} class.
 */
@Documented
@Inherited
@Retention(RetentionPolicy.RUNTIME)
@Target(ElementType.TYPE)
public @interface OperatorMetadata {
    /**
     * @return An alias name for the operator.
     */
    String alias();

    /**
     * @return The version of the operator.
     *         Defaults to the empty string (= not set).
     */
    String version() default "";

    /**
     * @return The author(s) of the operator.
     *         Defaults to the empty string (= not set).
     */
    String authors() default "";

    /**
     * @return The copyright notice for the operator code.
     *         Defaults to the empty string (= not set).
     */
    String copyright() default "";

    /**
     * @return A brief description of the operator's purpose.
     *         Defaults to the empty string (= not set).
     */
    String description() default "";

    /**
     * @return A Category to group the operator in.
     *         Defaults to the empty string (= not set).
     */
    String category() default "";
    
    /**
     * @return If {@code true}, this operator is considered for internal use onlyand thus
     *         may not be exposed in user interfaces.
     */
    boolean internal() default false;
}
