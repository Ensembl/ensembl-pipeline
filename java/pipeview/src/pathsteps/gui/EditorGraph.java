/*
 * @(#)EditorGraph.java	1.0.5 06/01/03
 *
 * Copyright (C) 2002 Gaudenz Alder
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
 */

package pathsteps.gui;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.event.MouseEvent;
import java.util.Hashtable;
import java.util.Map;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;
import javax.swing.event.CellEditorListener;

import org.jgraph.JGraph;
import org.jgraph.graph.CellView;
import org.jgraph.graph.GraphCellEditor;
import org.jgraph.graph.GraphConstants;
import org.jgraph.graph.GraphModel;
import org.jgraph.plaf.basic.BasicGraphUI;

/** 
* An example that demonstrates how to use a JDialog 
* as a CellEditor in JGraph. 
* 
* @version 1.1 23/12/02 
* @author Gaudenz Alder 
*/

public class EditorGraph extends JGraph {

	/** 
	* Constructs a EditorGraph with a sample model. 
	*/
	public EditorGraph() {
	}

	/** 
	* Constructs a EditorGraph for <code>model</code>. 
	*/
	public EditorGraph(GraphModel model) {
		super(model);
	}

	/** 
	* Override parent method with custom GraphUI. 
	*/
	public void updateUI() {
		setUI(new EditorGraphUI());
		invalidate();
	}

	public static void main(String[] args) {
		JFrame frame = new JFrame("EditorGraph");
		frame.getContentPane().add(new EditorGraph());
		frame.pack();
		//frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.setVisible(true);
	}

	/** 
	* Definition of the custom GraphUI. 
	*/
	public class EditorGraphUI extends BasicGraphUI {

		protected CellEditorListener cellEditorListener;

		protected JFrame editDialog = null;

		/** 
		* Create the dialog using the cell's editing component. 
		*/
		protected void createEditDialog(Object cell, MouseEvent event) {
			Dimension editorSize = editingComponent.getPreferredSize();
			editDialog = new JFrame("Edit " + graph.convertValueToString(cell));
			editDialog.setSize(editorSize.width, editorSize.height);
			editDialog.getContentPane().add(editingComponent);
			editingComponent.validate();
			editDialog.pack();
			editDialog.show();
		}

		/** 
		* Stops the editing session. If messageStop is true the editor 
		* is messaged with stopEditing, if messageCancel is true the 
		* editor is messaged with cancelEditing. If messageGraph is true 
		* the graphModel is messaged with valueForCellChanged. 
		*/
		protected void completeEditing(
			boolean messageStop,
			boolean messageCancel,
			boolean messageGraph) {
			if (stopEditingInCompleteEditing
				&& editingComponent != null
				&& editDialog != null) {
				Component oldComponent = editingComponent;
				Object oldCell = editingCell;
				GraphCellEditor oldEditor = cellEditor;
				Object newValue = oldEditor.getCellEditorValue();
				Rectangle editingBounds = graph.getCellBounds(editingCell);
				boolean requestFocus =
					(graph != null
						&& (graph.hasFocus() || editingComponent.hasFocus()));
				editingCell = null;
				editingComponent = null;
				if (messageStop)
					oldEditor.stopCellEditing();
				else if (messageCancel)
					oldEditor.cancelCellEditing();
				editDialog.dispose();
				if (requestFocus)
					graph.requestFocus();
				if (messageGraph) {
					Map map = GraphConstants.createMap();
					GraphConstants.setValue(map, newValue);
					Map nested = new Hashtable();
					nested.put(oldCell, map);
					graphLayoutCache.edit(nested, null, null, null);
				}
				updateSize();
				// Remove Editor Listener 
				if (oldEditor != null && cellEditorListener != null)
					oldEditor.removeCellEditorListener(cellEditorListener);
				cellEditor = null;
				editDialog = null;
			}
		}

		/** 
		* Will start editing for cell if there is a cellEditor and 
		* shouldSelectCell returns true.<p> 
		* This assumes that cell is valid and visible. 
		*/
		protected boolean startEditing(Object cell, MouseEvent event) {
			completeEditing();
			if (graph.isCellEditable(cell) && editDialog == null) {

				// Create Editing Component **** ***** 
				CellView tmp = graphLayoutCache.getMapping(cell, false);
				cellEditor = tmp.getEditor();
				editingComponent =
					cellEditor.getGraphCellEditorComponent(
						graph,
						cell,
						graph.isCellSelected(cell));
				if (cellEditor.isCellEditable(event)) {
					editingCell = cell;

					// Create Wrapper Dialog **** ***** 
					createEditDialog(cell, event);

					// Add Editor Listener 
					if (cellEditorListener == null)
						cellEditorListener = createCellEditorListener();
					if (cellEditor != null && cellEditorListener != null)
						cellEditor.addCellEditorListener(cellEditorListener);

					if (cellEditor.shouldSelectCell(event)) {
						stopEditingInCompleteEditing = false;
						try {
							graph.setSelectionCell(cell);
						} catch (Exception e) {
							System.err.println("Editing exception: " + e);
						}
						stopEditingInCompleteEditing = true;
					}

					if (event instanceof MouseEvent) {
						/* Find the component that will get forwarded all the 
						mouse events until mouseReleased. */
						Point componentPoint =
							SwingUtilities.convertPoint(
								graph,
								new Point(event.getX(), event.getY()),
								editingComponent);

						/* Create an instance of BasicTreeMouseListener to handle 
						passing the mouse/motion events to the necessary 
						component. */
						// We really want similiar behavior to getMouseEventTarget, 
						// but it is package private. 
						Component activeComponent =
							SwingUtilities.getDeepestComponentAt(
								editingComponent,
								componentPoint.x,
								componentPoint.y);
						if (activeComponent != null) {
							new MouseInputHandler(
								graph,
								activeComponent,
								event);
						}
					}
					return true;
				} else
					editingComponent = null;
			}
			return false;
		}

	}

}
