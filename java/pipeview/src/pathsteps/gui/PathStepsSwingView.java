package pathsteps.gui;
import pathsteps.common.*;
import pathsteps.model.*;
import java.awt.*;
import javax.swing.*;

import org.jgraph.*;
import org.jgraph.event.*;
import org.jgraph.graph.*;


/**
 * Swing implementation of the AView interface
**/
public class PathStepsSwingView implements AView
{
  
  private JFrame rootFrame;
  private JMenuBar menuBar;
  private JMenu fileMenu;
  private JMenuItem readNewPipelineDBMenuItem;
  private JMenuItem setNewLayoutConfigurationMenuItem;
  private JMenuItem saveGraphLayoutConfigurationMenuItem;
  
  private JMenu editMenu;
  private JMenuItem createMenuItem;
  private JMenuItem updateMenuItem;
  private JMenuItem deleteMenuItem;
  private JMenuItem connectMenuItem;
  
  private PathStepsPanel pathStepsPanel;
  private JTextField messageBar;
  private JButton saveButton;
  private JButton resetButton;
  private JButton refreshButton;
  private JButton layoutButton;
  private Application application;
  private ReadPipelineDBDialog readPipelineDBDialog;
  private LayoutConfigurationDialog layoutConfigurationDialog;
  private CreateNodeDialog createNodeDialog;
  private ALogger logger;
  
  public static String MAXIMUM_SCROLLPANE_WIDTH = "MAXIMUM_SCROLLPANE_WIDTH";
  public static String MAXIMUM_SCROLLPANE_HEIGHT = "MAXIMUM_SCROLLPANE_HEIGHT";
  
  public static String MINIMUM_SCROLLPANE_WIDTH = "MINIMUM_SCROLLPANE_WIDTH";
  public static String MINIMUM_SCROLLPANE_HEIGHT = "MINIMUM_SCROLLPANE_HEIGHT";
    
  public PathStepsSwingView(){}
  
  public void initialise(){
    int row = 0;

    if(getApplication() == null){
      throw new FatalAException("Application not set: Can't initialise view without Application");
    }
    
    setRootFrame(new JFrame("Pipe View"));
    setPathStepsPanel(new PathStepsPanel(this));
    setPathStepsPanel(getPathStepsPanel());
    setFileMenu(new JMenu("File"));
    setMenuBar(new JMenuBar());
    setReadNewPipelineDBMenuItem(new JMenuItem("Read New Pipeline DB"));
    setSetNewLayoutConfigurationMenuItem(new JMenuItem("Set New Layout Configuration"));
    setSaveGraphLayoutConfigurationMenuItem(new JMenuItem("Save Graph Layout Configuration"));

    setEditMenu(new JMenu("Edit"));
    setCreateMenuItem(new JMenuItem("Create Node"));
    setUpdateMenuItem(new JMenuItem("Update Node"));
    setDeleteMenuItem(new JMenuItem("Delete Node"));
    setConnectMenuItem(new JMenuItem("Connect Nodes"));
    
    connectToWindowCloseEventRouter(
      getRootFrame(),
      CLOSE_ROOT_FRAME_KEY
    );
    
    getRootFrame().setJMenuBar(getMenuBar());
    getMenuBar().add(getFileMenu());
    getFileMenu().add(getReadNewPipelineDBMenuItem());
    getFileMenu().add(getSetNewLayoutConfigurationMenuItem());
    getFileMenu().add(getSaveGraphLayoutConfigurationMenuItem());
    
    connectToActionEventRouter(
      getReadNewPipelineDBMenuItem(),
      AView.READ_NEW_PIPELINE_DB_KEY
    );
    
    connectToActionEventRouter(
      getSetNewLayoutConfigurationMenuItem(),
      AView.SET_NEW_LAYOUT_PREFERENCES_KEY
    );
    
    connectToActionEventRouter(
      getSaveGraphLayoutConfigurationMenuItem(),
      AView.SAVE_GRAPH_LAYOUT_CONFIGURATION_KEY
    );
    
    /* SWITCH ON EDITING SHORTLY...
    getMenuBar().add(getEditMenu());
    getEditMenu().add(getCreateMenuItem());
    getEditMenu().add(getUpdateMenuItem());
    getEditMenu().add(getDeleteMenuItem());
    getEditMenu().add(getConnectMenuItem());
    
    connectToActionEventRouter(
      getCreateMenuItem(),
      AView.CREATE_NODE_DIALOG_KEY
    );
    
    connectToActionEventRouter(
      getUpdateMenuItem(),
      AView.UPDATE_NODE_DIALOG_KEY
    );
    
    connectToActionEventRouter(
      getDeleteMenuItem(),
      AView.DELETE_NODE_DIALOG_KEY
    );
    
    connectToActionEventRouter(
      getConnectMenuItem(),
      AView.CONNECT_NODE_DIALOG_KEY
    );
    */
    
    setMessageBar(new JTextField(50));
    getMessageBar().setEditable(false);
    getMessageBar().setMinimumSize(new Dimension(600, 25));
    
    setSaveButton(new JButton("Save"));
    connectToActionEventRouter(
      getSaveButton(),
      SAVE_BUTTON_KEY
    );
    
    setResetButton(new JButton("Reset"));
    connectToActionEventRouter(
      getResetButton(),
      RESET_BUTTON_KEY
    );
    
    setRefreshButton(new JButton("Refresh"));
    connectToActionEventRouter(
      getRefreshButton(),
      REFRESH_BUTTON_KEY
    );
    
    setLayoutButton(new JButton("Layout"));
    connectToActionEventRouter(
      getLayoutButton(),
      LAYOUT_BUTTON_KEY
    );
    
    getRootFrame().getContentPane().setLayout(new GridBagLayout());
    
    GridBagConstraints constraint = makeConstraintAt(0, row, 4);
    constraint.fill = GridBagConstraints.BOTH;
    getPathStepsPanel().setMinimumSize(new Dimension(500, 500));
    getRootFrame().getContentPane().add(
      getPathStepsPanel(), 
      constraint
    );
    
    /*
    row++;

    getRootFrame().getContentPane().add(
      getSaveButton(), 
      makeConstraintAt(0,row,1)
    );
    
    getRootFrame().getContentPane().add(
      getResetButton(), 
      makeConstraintAt(1,row,1)
    );

     */
    
    row++;
    
    constraint = makeConstraintAt(0,row,2);
    constraint.fill = GridBagConstraints.BOTH;
    getRootFrame().getContentPane().add(
      getMessageBar(), 
      constraint
    );
    
    constraint = makeConstraintAt(2,row,1);
    constraint.anchor = GridBagConstraints.WEST;
    getRootFrame().getContentPane().add(
      getRefreshButton(), 
      constraint
    );

    constraint = makeConstraintAt(3,row,1);
    constraint.anchor = GridBagConstraints.WEST;
    getRootFrame().getContentPane().add(
      getLayoutButton(), 
      constraint
    );

  }

  public void show(){
    getRootFrame().pack();
    getRootFrame().show();
  }
  
  public void reshow(){
    getRootFrame().getContentPane().remove(getPathStepsPanel());
    
    GridBagConstraints constraint = makeConstraintAt(0, 0, 4);
    constraint.fill = GridBagConstraints.BOTH;
    getRootFrame().getContentPane().add(
      getPathStepsPanel(), 
      constraint
    );

    getRootFrame().invalidate();
    getRootFrame().validate();
    getRootFrame().repaint();
  }
  
  /**
   * Read the model and put its contents into the view.
  **/
  public void read(AModel aModel){
    PathStepsModel model = (PathStepsModel)aModel;
    if(getPathStepsPanel() != null){
      getPathStepsPanel().read(model);
    }
    
    if(getReadPipelineDBDialog() != null){
      getReadPipelineDBDialog().read(model);
    }
    
    if(getLayoutConfigurationDialog() != null){
      getLayoutConfigurationDialog().read(model);
    }
    
    if(getCreateNodeDialog() != null){
      getCreateNodeDialog().read(model);
    }
    
    showMessage((String)model.getRootElement().getProperty(PathStepsModel.MESSAGE));
    
    getRootFrame().invalidate();
    getRootFrame().validate();
    getRootFrame().repaint();
  }

  public void applyGraphLayout(){
    getPathStepsPanel().applyGraphLayout();
  }
  
  /**
   * update the model to reflect the view.
  **/
  public void update(AModel aModel) {
    PathStepsModel model = (PathStepsModel)aModel;

    if(getReadPipelineDBDialog() != null){
      getReadPipelineDBDialog().update(model);
    }
    
    if(getLayoutConfigurationDialog() != null){
      getLayoutConfigurationDialog().update(model);
    }
    
    if(getPathStepsPanel() != null){
      getPathStepsPanel().update(model);
    }
  }
  
  public Application getApplication(){
    return application;
  }
  
  public void setApplication(Application application){
    this.application = application;
  }
  
  public void showMessage(String message){
    getMessageBar().setText(message);
  }
  

  public void showFatalDialog(String message){
    JOptionPane.showMessageDialog(getRootFrame(), "Fatal Application Problem\n"+message);
  }
  
  public void shutDown(){
    getRootFrame().setVisible(false);
    getRootFrame().dispose();
    System.exit(0);
  }
  
  ////////
  ////////
  ////////
  public void openReadPipelineDBDialog(){
    setReadPipelineDBDialog(new ReadPipelineDBDialog(this));
    getReadPipelineDBDialog().pack();
    getReadPipelineDBDialog().show();
  }

  
  public void closeReadPipelineDBDialog() {
    getReadPipelineDBDialog().setVisible(false);
    getReadPipelineDBDialog().dispose();
    setReadPipelineDBDialog(null);
  }
  
  public boolean isReadPipelineDBDialogOpen() {
    return getReadPipelineDBDialog() != null;
  }
  
  public void requestFocusOnPipelineDBDialog() {
    if(!isReadPipelineDBDialogOpen()){
      throw new FatalAException("Attempt to request focus on ReadPipelineDBDialog when it's not open");
    }
    getReadPipelineDBDialog().requestFocus();
  }
  
  public void openLayoutConfigurationDialog(){
    setLayoutConfigurationDialog(new LayoutConfigurationDialog(this));
    getLayoutConfigurationDialog().pack();
    getLayoutConfigurationDialog().show();
  }

  
  public void closeLayoutConfigurationDialog(){
    getLayoutConfigurationDialog().setVisible(false);
    getLayoutConfigurationDialog().dispose();
    setLayoutConfigurationDialog(null);
  }
  
  public boolean isLayoutConfigurationDialogOpen() {
    return getLayoutConfigurationDialog() != null;
  }
  
  public void requestFocusOnLayoutConfigurationDialog() {
    if(!isLayoutConfigurationDialogOpen()){
      throw new FatalAException("Attempt to request focus on ReadPipelineDBDialog when it's not open");
    }
    getLayoutConfigurationDialog().requestFocus();
  }
  
  public ALogger getLogger() {
    return logger;
  }
  
  public void setLogger(ALogger logger) {
    this.logger = logger;
  }
  
  ////////
  ////////
  ////////
  JFrame getRootFrame(){
    return rootFrame;
  }
  
  private void setRootFrame(JFrame frame){
    rootFrame = frame;
  }
  
  private PathStepsPanel getPathStepsPanel(){
    return pathStepsPanel;
  }
  
  private void setPathStepsPanel(PathStepsPanel panel){
    pathStepsPanel = panel;
  }
  
  private JTextField getMessageBar(){
    return messageBar;
  }
  
  private void setMessageBar(JTextField field){
    messageBar = field;
  }
  
  private JButton getSaveButton(){
    return saveButton;
  }
  
  private void setSaveButton(JButton button){
    saveButton = button;
  }
  
  private JButton getResetButton(){
    return resetButton;
  }
  
  private void setResetButton(JButton button){
    resetButton = button;
  }

  private JButton getLayoutButton(){
    return layoutButton;
  }
  
  private void setLayoutButton(JButton button){
    layoutButton = button;
  }

  private JButton getRefreshButton(){
    return refreshButton;
  }
  
  private void setRefreshButton(JButton button){
    refreshButton = button;
  }
  
  private GridBagConstraints makeConstraintAt(
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
  
  private JMenuBar getMenuBar(){
    return menuBar;
  }

  private void setMenuBar(JMenuBar menuBar){
    this.menuBar = menuBar;
  }
  
  private JMenu getFileMenu(){
    return fileMenu;
  }

  private void setFileMenu(JMenu menu){
    fileMenu = menu;
  }
  
  private JMenuItem getReadNewPipelineDBMenuItem(){
    return readNewPipelineDBMenuItem;
  }
  
  private void setReadNewPipelineDBMenuItem(JMenuItem menuItem){
    readNewPipelineDBMenuItem = menuItem;
  }
  
  private void setSetNewLayoutConfigurationMenuItem(JMenuItem menuItem){
    setNewLayoutConfigurationMenuItem = menuItem;
  }
  
  private JMenuItem getSetNewLayoutConfigurationMenuItem(){
    return setNewLayoutConfigurationMenuItem;
  }

  private JMenu getEditMenu(){
    return editMenu;
  }
  
  private void setEditMenu(JMenu newValue){
    editMenu = newValue;
  }
  
  private JMenuItem getCreateMenuItem(){
    return createMenuItem;
  }
  
  private void setCreateMenuItem(JMenuItem newValue){
    createMenuItem = newValue;
  }
  
  private JMenuItem getUpdateMenuItem(){
    return updateMenuItem;
  }
  
  private void setUpdateMenuItem(JMenuItem newValue){
    updateMenuItem = newValue;
  }
  
  private JMenuItem getDeleteMenuItem(){
    return deleteMenuItem;
  }
  
  private void setDeleteMenuItem(JMenuItem newValue){
    deleteMenuItem = newValue;
  }
  
  private JMenuItem getConnectMenuItem(){
    return connectMenuItem;
  }  
  
  private void setConnectMenuItem(JMenuItem newValue){
    connectMenuItem = newValue;
  }
  
  private JMenuItem getSaveGraphLayoutConfigurationMenuItem(){
    return saveGraphLayoutConfigurationMenuItem;
  }
  
  private void setSaveGraphLayoutConfigurationMenuItem(JMenuItem item){
    saveGraphLayoutConfigurationMenuItem = item;
  }
  
  private ReadPipelineDBDialog getReadPipelineDBDialog(){
    return readPipelineDBDialog;
  }

  private void setReadPipelineDBDialog(ReadPipelineDBDialog dialog){
    readPipelineDBDialog = dialog;
  }

  private LayoutConfigurationDialog getLayoutConfigurationDialog(){
    return layoutConfigurationDialog ;
  }
  
  private void setLayoutConfigurationDialog(LayoutConfigurationDialog dialog){
    layoutConfigurationDialog = dialog;
  }
  
  private CreateNodeDialog getCreateNodeDialog(){
    return createNodeDialog;
  }
  
  private void setCreateNodeDialog(CreateNodeDialog dialog){
    createNodeDialog = dialog;
  }
  
  void connectToActionEventRouter(
    JButton button,
    String key
  ){
    button.addActionListener(
      new ActionEventRouter(
        getApplication(),
        key
      )
    );
  }
  
  void connectToActionEventRouter(
    JMenuItem item,
    String key
  ){
    item.addActionListener(
      new ActionEventRouter(
        getApplication(),
        key
      )
    );
  }
  
  void connectToWindowCloseEventRouter(
    java.awt.Window frame,
    String key
  ){
    frame.addWindowListener(
      new WindowCloseEventRouter(
        getApplication(),
        key
      )
    );
  }
  
  void connectToActionEventRouter(
    JComboBox item,
    String key
  ){
    item.addActionListener(
      new ActionEventRouter(
        getApplication(),
        key
      )
    );
  }

  void connectToKeyEventRouter(
    JTextField item,
    String key
  ){
    item.addKeyListener(
      new KeyEventRouter(
        getApplication(),
        key
      )
    );
  }
  
  void connectToMouseLeftClickRouter(
    JGraph graph,
    String key
  ){
    graph.addMouseListener(
      new MouseLeftClickEventRouter(
        getApplication(),
        key
      )
    );
  }
  
  public void clear(){
    getPathStepsPanel().clear();
  }
  
  public Object getSelectedGraphCellObject(){
    if(getPathStepsPanel()==null){
      throw new FatalAException("Attempt made to find selected graph node when panel not initialised");
    }
    
    return getPathStepsPanel().getSelectedGraphCellObject();
  }
  
  public void closeCreateNodeDialog() {
    getCreateNodeDialog().setVisible(false);
    getCreateNodeDialog().dispose();
    getCreateNodeDialog();
  }
  
  public boolean isCreateNodeDialogOpen() {
    return getCreateNodeDialog() != null;
  }
  
  public void openCreateNodeDialog() {
    setCreateNodeDialog(new CreateNodeDialog(this));
    getCreateNodeDialog().pack();
    getCreateNodeDialog().show();
  }
  
  public void requestFocusOnCreateNodeDialog() {
      if(!isCreateNodeDialogOpen()){
      throw new FatalAException("Attempt to request focus on create node dialog when it's not open");
    }
    getCreateNodeDialog().requestFocus();
  }
  
  public String getNameOfSelectedNode() {
    DefaultGraphCell cell;
    Object value = null;
    Object selected = getSelectedGraphCellObject();
    if(
      selected instanceof DefaultEdge
    ){
      if(getLogger().isLoggingMedium()){
        getLogger().logMedium("Edge selected - doing nothing");
      }
      
      return null;
      
    }else if(
      selected instanceof DefaultGraphCell
    ){
      cell = (DefaultGraphCell)selected;
      value = GraphConstants.getValue(cell.getAttributes());
      
      if(!(value instanceof String)){
        throw new FatalAException("Expected a String for cell's value object: found "+value);
      }

      return (String)value;
      
    }else{
      throw new FatalAException("Something was selected: "+selected+" and I don't know what it is");
    }
    

  }

  public boolean isLoggingLow(){
    return getLogger().isLoggingLow();
  }
  
  public void logLow(String message){
    getLogger().logLow(message);
  }
  
  public void logLow(String message, Throwable exception) {
    getLogger().logLow(message, exception);
  }
  
  public boolean isLoggingMedium() {
    return getLogger().isLoggingMedium();
  }
  
  public void logMedium(String message) {
    getLogger().logMedium(message);
  }
  
  public void logMedium(String message, Throwable exception) {
    getLogger().logMedium(message, exception);
  }
  
  public boolean isLoggingHigh(){
    return getLogger().isLoggingHigh();
  }
  
  public void logHigh(String message) {
    getLogger().logHigh(message);
  }
  
  public void logHigh(String message, Throwable exception) {
    getLogger().logHigh(message, exception);
  }
}