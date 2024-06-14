/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2014-2018 The Contributors of the realKD Project
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
package de.unibonn.realkd.patterns.correlated;

import static de.unibonn.realkd.common.base.Identifier.id;
import static de.unibonn.realkd.data.Populations.population;
import static de.unibonn.realkd.data.table.DataTables.table;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;

import de.unibonn.realkd.data.Population;
import de.unibonn.realkd.data.table.DataTable;
import de.unibonn.realkd.data.table.attribute.Attributes;
import de.unibonn.realkd.data.table.attribute.CategoricAttribute;

/**
 * @author Panagiotis Mandros
 *
 */
public class CorrelatedPatternTest {

	private Population population = population(id("test_population"), 4);

	private CategoricAttribute<String> a = Attributes.categoricalAttribute("A", "",
			ImmutableList.of("a", "a", "b", "b"));

	private CategoricAttribute<String> b = Attributes.categoricalAttribute("B", "",
			ImmutableList.of("a", "b", "a", "b"));

	private CategoricAttribute<String> c = Attributes.categoricalAttribute("C", "",
			ImmutableList.of("a", "a", "b", "b"));

	private DataTable table = table(id("table"), "table", "", population, ImmutableList.of(a, b, c));

	@Test
	public void constructionOfAttributeSetRelationTest() {
		AttributeSetRelation relation1 = CorrelationPatterns.attributeSetRelation(table, ImmutableSet.of(a, b, c));
		assertNotNull(relation1);
		assertEquals("Set must be {A,B,C}.", relation1.attributeSet(), ImmutableSet.of(a, b, c));

		AttributeSetRelation relation2 = CorrelationPatterns.attributeSetRelation(table, ImmutableSet.of(c, b, a));
		assertEquals("Set must be {A,B,C}.", relation2.attributeSet(), ImmutableSet.of(a, b, c));

	}

	@Test
	public void constructionOfCorrelationPatternTest() {
		AttributeSetRelation relation = CorrelationPatterns.attributeSetRelation(table, ImmutableSet.of(a, b, c));
		CorrelationPattern pattern = CorrelationPatterns.correlationPattern(relation);
		assertNotNull(pattern);
		assertEquals(pattern.descriptor(), relation);
		assertTrue(pattern.hasMeasure(NormalizedTotalCorrelation.NORMALIZED_TOTAL_CORRELATION));
	}
}
