/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2014-15 The Contributors of the realKD Project
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
package de.unibonn.realkd.patterns;

import de.unibonn.realkd.common.workspace.HasSerialForm;

/**
 * <p>
 * Syntactic expression that refers to some abstract concept within one or more
 * data artifacts---usually based on the meta-data of these artifacts. Examples:
 * </p>
 * <ul>
 * <li>"The propositions <em>x</em>, <em>y</em>, <em>z</em>"
 * <li>"The values of attributes <em>a</em>, <em>b</em> on the subset of data
 * described by the propositions <em>x</em>, <em>y</em>, <em>z</em>"
 * <li>"The rule that propositions <em>x</em>, <em>y</em> imply proposition
 * <em>z</em>"
 * </ul>
 * 
 * @author Mario Boley
 * 
 * @since 0.1.0
 * 
 * @version 0.1.2.1
 *
 */
public interface PatternDescriptor extends
		HasSerialForm<PatternDescriptor> {

}
