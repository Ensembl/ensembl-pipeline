package pathsteps.gui;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.util.*;

import pathsteps.common.*;
import pathsteps.model.*;

public class 
  DataSourceConfigurationPanel 
extends 
  JPanel
{
  private JLabel hostLabel = new JLabel("Host");
  private JLabel portLabel = new JLabel("Port");
  private JLabel userLabel = new JLabel("User");
  private JLabel passwordLabel = new JLabel("Password");
  private JLabel ensemblDatabaseLabel = new JLabel("Ensembl Database Name");
  
  private JTextField hostTextField = new JTextField(40);
  private JTextField portTextField = new JTextField(6);
  private JTextField userTextField = new JTextField(40);
  private JPasswordField passwordTextField = new JPasswordField(40);
  private JComboBox ensemblDatabaseDropdown = new JComboBox();
  
  private JButton findButton = new JButton("Find DBs...");
  private JButton readButton = new JButton("Read DB");
  private JButton cancelButton = new JButton("Cancel");
  
  private PathStepsSwingView view;
  
  public DataSourceConfigurationPanel(PathStepsSwingView view){
    this.view = view;
    buildGUI();
    //setDefaults();
  }
  
  private void setDefaults(){
    getHostTextField().setText("ecs1b");
    getPortTextField().setText("3306");
    getUserTextField().setText("ensro");
    getPasswordTextField().setText(null);
    getEnsemblDatabaseDropdown().setSelectedItem("vivek_pseudo_test");
  }

  private void buildGUI(){
    int row = 0;
    
    setLayout(new GridBagLayout());

    add(hostLabel,makeConstraintAt(0,row,1));
    add(hostTextField,makeConstraintAt(1,row,3));
    getView().connectToKeyEventRouter(
      hostTextField,
      getView().DATABASE_DATA_CHANGE_KEY
    );
    
    row++;
    add(portLabel,makeConstraintAt(0,row,1));
    add(portTextField,makeConstraintAt(1,row,3));
    getView().connectToKeyEventRouter(
      portTextField,
      getView().DATABASE_DATA_CHANGE_KEY
    );
    
    row++;
    add(userLabel,makeConstraintAt(0,row,1));
    add(userTextField,makeConstraintAt(1,row,3));
    getView().connectToKeyEventRouter(
      userTextField,
      getView().DATABASE_DATA_CHANGE_KEY
    );
    
    row++;
    add(passwordLabel,makeConstraintAt(0,row,1));
    add(passwordTextField,makeConstraintAt(1,row,3));
    getView().connectToKeyEventRouter(
      passwordTextField,
      getView().DATABASE_DATA_CHANGE_KEY
    );
    
    row++;
    add(ensemblDatabaseLabel,makeConstraintAt(0,row,1));

    //I want the ensembl db-dropdown to be 2/3 the size of the other textfields.
    Dimension hostTextFieldDimension = hostTextField.getPreferredSize();
    int newX = new Double((new Integer(hostTextFieldDimension.width).doubleValue())*0.66).intValue();
    Dimension newDimension = new Dimension(newX, hostTextFieldDimension.height);
    ensemblDatabaseDropdown.setPreferredSize(newDimension);
    add(ensemblDatabaseDropdown,makeConstraintAt(1,row,1));
    
    getView().connectToActionEventRouter(
      findButton,
      getView().FIND_ALL_DATABASES_KEY
    );
    
    add(findButton,makeConstraintAt(2,row,1));
    
    row++;
    getView().connectToActionEventRouter(
      readButton,
      getView().READ_SPECIFIED_DB_KEY
    );
    
    add(readButton,makeConstraintAt(0,row,1));
    
    getView().connectToActionEventRouter(
      cancelButton,
      getView().CANCEL_READ_KEY
    );
    
    add(cancelButton,makeConstraintAt(1,row,1));
  }
  
  public void read(PathStepsModel model){
    if(getView().getLogger().isLoggingMedium()){
      getView().getLogger().logMedium("Starting data source configuration panel read");
    }
    String host = null;
    String port = null;
    String user = null;
    String password = null;
    Collection databases = null;
    String selectedDatabase = null;
    ModelElement dialogModel = model.getRootElement().getChildElement(model.READ_DB_DIALOG);
    Collection list;
    
    if(dialogModel == null){
      throw new FatalAException("ReadDBDialog expecting a model called: "+model.READ_DB_DIALOG+" but none was found");
    }
    
    host = (String)dialogModel.getProperty(model.READ_DB_DIALOG_HOST);
    port = (String)dialogModel.getProperty(model.READ_DB_DIALOG_PORT);
    user = (String)dialogModel.getProperty(model.READ_DB_DIALOG_USER);
    password = (String)dialogModel.getProperty(model.READ_DB_DIALOG_PASSWORD);
    selectedDatabase = (String)dialogModel.getProperty(model.READ_DB_DIALOG_PIPELINE_DB);
    
    list = 
      (ArrayList)dialogModel
        .getChildElement(model.READ_DB_DIALOG_PIPELINE_DB_LIST)
        .getProperty(model.READ_DB_DIALOG_PIPELINE_DB_LIST);
    
    if(getView().getLogger().isLoggingMedium()){
      getView().getLogger().logMedium("list: " + list);
    }
    
    getHostTextField().setText(host);
    getPortTextField().setText(port);
    getUserTextField().setText(user);
    getPasswordTextField().setText(password);
    getEnsemblDatabaseDropdown().setModel(new DefaultComboBoxModel(new Vector(list)));
    getEnsemblDatabaseDropdown().setSelectedItem(selectedDatabase);
    if(getView().getLogger().isLoggingMedium()){
      getView().getLogger().logMedium("Finished data source configuration panel read");
    }
  }
  
  public void update(PathStepsModel model){
    String host = getHostTextField().getText();;
    String port = getPortTextField().getText();
    String user = getUserTextField().getText();
    String password = new String(getPasswordTextField().getPassword());
    String selectedDatabase = (String)getEnsemblDatabaseDropdown().getSelectedItem();

    ModelElement dialogModel = model.getRootElement().getChildElement(model.READ_DB_DIALOG);

    if(dialogModel == null){
      throw new FatalAException("ReadDBDialog expecting a model called: "+model.READ_DB_DIALOG+" but none was found");
    }
    
    dialogModel.addProperty(model.READ_DB_DIALOG_HOST,host);
    dialogModel.addProperty(model.READ_DB_DIALOG_PORT, port);
    dialogModel.addProperty(model.READ_DB_DIALOG_USER, user);
    dialogModel.addProperty(model.READ_DB_DIALOG_PASSWORD, password);
    dialogModel.addProperty(model.READ_DB_DIALOG_PIPELINE_DB, selectedDatabase);
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

  private JTextField getHostTextField(){
    return hostTextField;
  }
  
  private JTextField getPortTextField(){
    return portTextField;
  }
  
  private JTextField getUserTextField(){
    return userTextField;
  }
  
  private JPasswordField getPasswordTextField(){
    return passwordTextField;
  }
  
  public JComboBox getEnsemblDatabaseDropdown(){
    return ensemblDatabaseDropdown;
  }
  
  public String getSelectedEnsemblDatabase(){
    return (String)getEnsemblDatabaseDropdown().getSelectedItem();
  }
  
  public void setSelectedEnsemblDatabase(String database){
    getEnsemblDatabaseDropdown().setSelectedItem(database);
  }
  
  public JButton getFindButton(){
    return findButton;
  }
  
  public JButton getReadButton(){
    return readButton;
  }
  
  public JButton getCancelButton(){
    return cancelButton;
  }
  
  private PathStepsSwingView getView(){
    return view;
  }
}

