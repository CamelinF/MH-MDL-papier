/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2014-2017 The Contributors of the realKD Project
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
package de.unibonn.realkd.data.table;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import de.unibonn.realkd.common.base.Identifier;
import de.unibonn.realkd.data.xarf.XarfImport;

/**
 * Tests correct parsing of octet binary dataset from xarf file.
 * 
 * @author Mario Boley
 * 
 * @since 0.5.0
 * 
 * @version 0.5.0
 *
 */
public class BinariesXarfInputTest {

	public static final String FILENAME = "src/main/resources/data/binaries/octet_binaries_2.1.1.xarf";

	private final static XarfImport builder = XarfImport.xarfImport()
			.dataFilename(FILENAME);
	private final static DataTable dataTable = builder.get();

	@Test
	public void builtNonNullTable() {
		assertNotNull(dataTable);
	}

	@Test
	public void testNumberOfAttributes() {
		assertEquals(56, dataTable.attributes().size());
	}

	@Test
	public void testGroupParsed() {
		assertEquals(1, dataTable.attributeGroups().size());
	}

	@Test
	public void testNamesParsed() {
		assertEquals("BeO", dataTable.population().objectName(4));
	}

	@Test
	public void testTableNameParsed() {
		assertEquals("Octet Binaries", dataTable.caption());
	}
	
	@Test
	public void testAttributeIdentifiersAccess() {
		assertTrue(dataTable.attribute(Identifier.identifier("lasso_1")).isPresent());
		assertFalse(dataTable.attribute(Identifier.identifier("no_such_identifier")).isPresent());
	}

}