package pathsteps.gui;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.util.*;

import pathsteps.common.*;
import pathsteps.model.*;

public class LayoutConfigurationPanel extends JPanel{
  private JLabel iteratesLabel = new JLabel("Iterates");
  private JLabel springLengthLabel = new JLabel("Spring Length");
  private JLabel repulsionMultiplierLabel = new JLabel("Repulsion Multiplier");
  private JLabel horizontalSpacingLabel = new JLabel("Horizontal Spacing");
  private JLabel verticalSpacingLabel = new JLabel("Vertical Spacing");
  private JLabel movementLimitLabel = new JLabel("Movement Limit");
  private JLabel gravityLabel = new JLabel("Gravity");
  
  private JTextField iteratesTextField = new JTextField(10);
  private JTextField springLengthTextField = new JTextField(10);
  private JTextField repulsionMultiplierTextField = new JTextField(10);
  private JTextField horizontalSpacingTextField = new JTextField(10);
  private JTextField verticalSpacingTextField = new JTextField(10);
  private JTextField movementLimitTextField = new JTextField(10);
  private JTextField gravityTextField = new JTextField(10);
  private JCheckBox fixNodesCheckBox = new JCheckBox("Fix nodes");
  private JCheckBox showJobDetailCheckBox = new JCheckBox("Show job detail");
  
  private JButton setButton = new JButton("Set");
  private JButton cancelButton = new JButton("Cancel");
  
  private PathStepsSwingView view;
  
  public LayoutConfigurationPanel(PathStepsSwingView view){
    this.view = view;
    buildGUI();
    //setDefaults();
  }
  
  private void buildGUI(){
    int row = 0;
    
    setLayout(new GridBagLayout());

    add(iteratesLabel,makeConstraintAt(0,row,1));
    add(getIteratesTextField(),makeConstraintAt(1,row,1));
    row++;
    
    add(springLengthLabel,makeConstraintAt(0,row,1));
    add(getSpringLengthTextField(),makeConstraintAt(1,row,2));
    row++;
    
    add(repulsionMultiplierLabel,makeConstraintAt(0,row,1));
    add(getRepulsionMultiplierTextField(),makeConstraintAt(1,row,2));
    row++;
    
    add(horizontalSpacingLabel,makeConstraintAt(0,row,1));
    add(getHorizontalSpacingTextField(),makeConstraintAt(1,row,2));
    row++;

    add(verticalSpacingLabel,makeConstraintAt(0,row,1));
    add(getVerticalSpacingTextField(),makeConstraintAt(1,row,2));
    row++;

    add(movementLimitLabel,makeConstraintAt(0,row,1));
    add(getMovementLimitTextField(),makeConstraintAt(1,row,2));
    row++;

    add(gravityLabel,makeConstraintAt(0,row,1));
    add(getGravityTextField(),makeConstraintAt(1,row,2));
    row++;

    add(getFixNodesCheckBox(),makeConstraintAt(0,row,2));
    row++;

    add(getShowJobDetailCheckBox(),makeConstraintAt(0,row,2));
    row++;

    getView().connectToActionEventRouter(
      getSetButton(),
      AView.SET_LAYOUT_PREFERENCES_KEY
    );
    
    add(getSetButton(),makeConstraintAt(0,row,1));
    
    getView().connectToActionEventRouter(
      getCancelButton(),
      getView().CANCEL_LAYOUT_PREFERENCES_KEY
    );
    
    add(getCancelButton(),makeConstraintAt(1,row,1));
  }
  
  public void read(PathStepsModel model){
    if(getView().getLogger().isLoggingMedium()){
      getView().getLogger().logMedium("Starting layout configuration panel read");
    }
    
    String iterates;
    String springLength;
    String repulsionMultiplier;
    String horizontalSpacing;
    String verticalSpacing;
    String movementLimit;
    String gravity;
    String fixNodes;
    String showJobDetail;

    ModelElement dialogModel = model.getRootElement().getChildElement(model.LAYOUT_DIALOG);

    if(dialogModel == null){
      throw new FatalAException(
        "LayoutConfigurationDialog expecting a model called: "+
        model.LAYOUT_DIALOG+
        " but none was found"
      );
    }

    iterates = (String)dialogModel.getProperty(model.LAYOUT_DIALOG_ITERATES);
    springLength = (String)dialogModel.getProperty(model.LAYOUT_DIALOG_SPRING_NATURAL_LENGTH);
    repulsionMultiplier = (String)dialogModel.getProperty(model.LAYOUT_DIALOG_REPULSION_MULTIPLIER);
    horizontalSpacing = (String)dialogModel.getProperty(model.LAYOUT_DIALOG_HORIZONTAL_SPACING);
    verticalSpacing = (String)dialogModel.getProperty(model.LAYOUT_DIALOG_VERTICAL_SPACING);
    movementLimit = (String)dialogModel.getProperty(model.LAYOUT_DIALOG_MOVEMENT_LIMIT);
    gravity = (String)dialogModel.getProperty(model.LAYOUT_DIALOG_GRAVITY);
    fixNodes = (String)dialogModel.getProperty(model.LAYOUT_DIALOG_FIX_ROOTS);
    showJobDetail = (String)dialogModel.getProperty(model.LAYOUT_DIALOG_SHOW_JOB_DETAIL);

    getIteratesTextField().setText(iterates);
    getSpringLengthTextField().setText(springLength);
    getRepulsionMultiplierTextField().setText(repulsionMultiplier);
    getHorizontalSpacingTextField().setText(horizontalSpacing);
    getVerticalSpacingTextField().setText(verticalSpacing);
    getMovementLimitTextField().setText(movementLimit);
    getGravityTextField().setText(gravity);
    getFixNodesCheckBox().setSelected(Boolean.valueOf(fixNodes).booleanValue());
    getShowJobDetailCheckBox().setSelected(Boolean.valueOf(showJobDetail).booleanValue());
    
    if(getView().getLogger().isLoggingMedium()){
      getView().getLogger().logMedium("Finished data source configuration panel read");
    }
  }
  
  public void update(PathStepsModel model){
    String iterates = getIteratesTextField().getText();;
    String springLength = getSpringLengthTextField().getText();
    String repulsionMultiplier = getRepulsionMultiplierTextField().getText();
    String horizontalSpacing = getHorizontalSpacingTextField().getText();
    String verticalSpacing = getVerticalSpacingTextField().getText();
    String movementLimit = getMovementLimitTextField().getText();
    String gravity = getGravityTextField().getText();
    String fixNodes = Boolean.valueOf(getFixNodesCheckBox().isSelected()).toString();
    String showJobDetail = Boolean.valueOf(getShowJobDetailCheckBox().isSelected()).toString();

    ModelElement dialogModel = model.getRootElement().getChildElement(model.LAYOUT_DIALOG);

    if(dialogModel == null){
      throw new FatalAException("ReadDBDialog expecting a model called: "+model.READ_DB_DIALOG+" but none was found");
    }

    dialogModel.addProperty(model.LAYOUT_DIALOG_ITERATES,iterates);
    dialogModel.addProperty(model.LAYOUT_DIALOG_SPRING_NATURAL_LENGTH, springLength);
    dialogModel.addProperty(model.LAYOUT_DIALOG_REPULSION_MULTIPLIER, repulsionMultiplier);
    dialogModel.addProperty(model.LAYOUT_DIALOG_HORIZONTAL_SPACING, horizontalSpacing);
    dialogModel.addProperty(model.LAYOUT_DIALOG_VERTICAL_SPACING, verticalSpacing);
    dialogModel.addProperty(model.LAYOUT_DIALOG_MOVEMENT_LIMIT, movementLimit);
    dialogModel.addProperty(model.LAYOUT_DIALOG_GRAVITY, gravity);
    dialogModel.addProperty(model.LAYOUT_DIALOG_FIX_ROOTS, fixNodes);
    dialogModel.addProperty(model.LAYOUT_DIALOG_SHOW_JOB_DETAIL, showJobDetail);
  }

  protected GridBagConstraints makeConstraintAt(
    int x,
    int y,
    int width
  ) {
    GridBagConstraints gbc = new GridBagConstraints();
    gbc.gridx = x;
    gbc.gridy = y;
    gbc.gridwidth = width;
    gbc.gridheight = 1;
    gbc.weightx = 0.0;
    gbc.weighty = 0.0;
    gbc.anchor = GridBagConstraints.WEST;
    gbc.fill = GridBagConstraints.NONE;
    gbc.insets = new Insets(0, 0, 0, 0);
    return gbc;
  }//end makeConstraintAt

  private JTextField getIteratesTextField(){
    return iteratesTextField;
  }
  
  private JTextField getSpringLengthTextField(){
    return springLengthTextField;
  }
  
  private JTextField getRepulsionMultiplierTextField(){
    return repulsionMultiplierTextField;
  }
  
  private JTextField getHorizontalSpacingTextField(){
    return horizontalSpacingTextField;
  }
  
  private JTextField getVerticalSpacingTextField(){
    return verticalSpacingTextField;
  }
  
  private JTextField getMovementLimitTextField(){
    return movementLimitTextField;
  }
  
  private JTextField getGravityTextField(){
    return gravityTextField;
  }
  
  public JButton getSetButton(){
    return setButton;
  }
  
  public JButton getCancelButton(){
    return cancelButton;
  }
  
  private PathStepsSwingView getView(){
    return view;
  }

  private JCheckBox getShowJobDetailCheckBox(){
    return showJobDetailCheckBox;
  }
  
  private JCheckBox getFixNodesCheckBox(){
    return fixNodesCheckBox;
  }
}

