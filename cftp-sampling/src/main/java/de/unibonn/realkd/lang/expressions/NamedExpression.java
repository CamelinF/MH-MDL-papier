/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2014-16 The Contributors of the realKD Project
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 */
package de.unibonn.realkd.lang.expressions;

import de.unibonn.realkd.common.workspace.Workspace;

/**
 * @author boley
 *
 */
public class NamedExpression<T> implements Expression<T> {

	private final String name;

	private final Expression<T> expression;

	private T value = null;

	public NamedExpression(String name, Expression<T> expression) {
		this.name = name;
		this.expression=expression;
	}

	@Override
	public T evaluate(Workspace workspace) throws InterpretationException {
		return value == null ? value = expression.evaluate(workspace) : value;
	}

	@Override
	public Class<?> resultType() {
		return expression.resultType();
	}
	
	public String name() {
		return name;
	}

}
