package pathsteps.gui;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.util.*;

import pathsteps.common.*;
import pathsteps.model.*;

public class CreateNodePanel extends JPanel{
  private JLabel nameLabel = new JLabel("Name");
  private JLabel typeLabel = new JLabel("Type");
  private JLabel dbLabel = new JLabel("Db");
  private JLabel dbFileLabel = new JLabel("DbFile");
  private JLabel programLabel = new JLabel("Program");
  private JLabel programVersionLabel = new JLabel("ProgramVersion");
  private JLabel programFileLabel = new JLabel("ProgramFile");
  private JLabel parametersLabel = new JLabel("Parameters");
  private JLabel moduleLabel = new JLabel("Module");
  private JLabel gffSourceLabel = new JLabel("GFF Source");
  private JLabel gffFeatureLabel = new JLabel("GFF Feature");

  private JTextField nameTextField = new JTextField(10);
  private JTextField dbTextField = new JTextField(10);
  private JTextField dbFileTextField = new JTextField(10);
  private JTextField programTextField = new JTextField(10);
  private JTextField programVersionTextField = new JTextField(10);
  private JTextField programFileTextField = new JTextField(10);
  private JTextField parametersTextField = new JTextField(10);
  private JTextField moduleTextField = new JTextField(10);
  private JTextField gffSourceTextField = new JTextField(10);
  private JTextField gffFeatureTextField = new JTextField(10);
  private JComboBox typeDropdown = new JComboBox();

  private JButton createButton = new JButton("Create");
  private JButton cancelButton = new JButton("Cancel");
  
  private PathStepsSwingView view;
  
  public CreateNodePanel(PathStepsSwingView view){
    this.view = view;
    buildGUI();
  }
  
  private void buildGUI(){
    int row = 0;
    
    setLayout(new GridBagLayout());

    add(nameLabel,makeConstraintAt(0,row,1));
    add(getNameTextField(),makeConstraintAt(1,row,1));
    row++;
    
    add(typeLabel,makeConstraintAt(0,row,1));
    add(getTypeDropdown(),makeConstraintAt(1,row,2));
    row++;
    
    add(dbFileLabel,makeConstraintAt(0,row,1));
    add(getDbFileTextField(),makeConstraintAt(1,row,1));
    row++;
    
    add(programLabel,makeConstraintAt(0,row,1));
    add(getProgramTextField(),makeConstraintAt(1,row,1));
    row++;
    
    add(programVersionLabel,makeConstraintAt(0,row,1));
    add(getProgramVersionTextField(),makeConstraintAt(1,row,1));
    row++;
    
    add(programFileLabel,makeConstraintAt(0,row,1));
    add(getProgramFileTextField(),makeConstraintAt(1,row,1));
    row++;
    
    add(parametersLabel,makeConstraintAt(0,row,1));
    add(getParametersTextField(),makeConstraintAt(1,row,1));
    row++;
    
    add(moduleLabel,makeConstraintAt(0,row,1));
    add(getModuleTextField(),makeConstraintAt(1,row,1));
    row++;
    
    add(gffSourceLabel,makeConstraintAt(0,row,1));
    add(getGffSourceTextField(),makeConstraintAt(1,row,1));
    row++;
    
    add(gffFeatureLabel,makeConstraintAt(0,row,1));
    add(getGffFeatureTextField(),makeConstraintAt(1,row,1));
    row++;
    
    getView().connectToActionEventRouter(
      getCreateButton(),
      AView.CONFIRM_CREATE_NODE_KEY
    );
    
    add(getCreateButton(),makeConstraintAt(0,row,1));
    
    getView().connectToActionEventRouter(
      getCancelButton(),
      getView().CANCEL_NODE_DIALOG_KEY
    );
    
    add(getCancelButton(),makeConstraintAt(1,row,1));
  }
  
  public void read(PathStepsModel model){
    if(getView().getLogger().isLoggingMedium()){
      getView().getLogger().logMedium("Starting layout configuration panel read");
    }
    
    ModelElement dialogModel = model.getRootElement().getChildElement(model.CREATE_NODE_DIALOG);
    String name;
    java.util.List allowedTypes;

    if(dialogModel == null){
      throw new FatalAException(
        "CreateNodeDialog expecting a model called: "+
        model.CREATE_NODE_DIALOG+
        " but none was found"
      );
    }

    name = (String)dialogModel.getProperty(model.CREATE_NODE_DIALOG_NAME);

    getNameTextField().setText(name);
    
    allowedTypes = (java.util.List)dialogModel.getProperty(model.CREATE_NODE_DIALOG_INPUT_ID_TYPE_LIST);

    if(allowedTypes == null){
      throw new FatalAException("Could not find 'allowedTypes' List in CreateNodePanel");
    }
    
    getTypeDropdown().setModel(
      new DefaultComboBoxModel(
        new Vector(allowedTypes)
      )
    );
    
    if(getView().getLogger().isLoggingMedium()){
      getView().getLogger().logMedium("Finished create node panel read");
    }
  }
  
  public void update(PathStepsModel model){
    String name = getNameTextField().getText();;

    ModelElement dialogModel = model.getRootElement().getChildElement(model.CREATE_NODE_DIALOG);

    if(dialogModel == null){
      throw new FatalAException("CreateNodeDialog expecting a model called: "+model.CREATE_NODE_DIALOG+" but none was found");
    }

    dialogModel.addProperty(model.CREATE_NODE_DIALOG,name);
    dialogModel.addProperty(model.CREATE_NODE_DIALOG_INPUT_ID_TYPE, (String)getTypeDropdown().getSelectedItem());
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

  private JTextField getNameTextField(){
    return nameTextField;
  }
  
  private JTextField getDbFileTextField(){
    return dbFileTextField;
  }
  
  private JTextField getProgramVersionTextField(){
    return programVersionTextField;
  }
  
  private JTextField getProgramTextField(){
    return programTextField;
  }
  
  private JTextField getProgramFileTextField(){
    return programFileTextField;
  }
  
  private JTextField getParametersTextField(){
    return parametersTextField;
  }
  
  private JTextField getModuleTextField(){
    return moduleTextField;
  }
  
  private JTextField getGffSourceTextField(){
    return gffSourceTextField;
  }
  
  private JTextField getGffFeatureTextField(){
    return gffFeatureTextField;
  }
    
  private JComboBox getTypeDropdown(){
    return typeDropdown;
  }
  
  public JButton getCreateButton(){
    return createButton;
  }
  
  public JButton getCancelButton(){
    return cancelButton;
  }
  
  private PathStepsSwingView getView(){
    return view;
  }
}

