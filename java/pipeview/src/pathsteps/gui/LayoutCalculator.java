package pathsteps.gui;
import java.awt.*;
import java.util.*;

public class LayoutCalculator {
  private int nnodes;
  private int nedges;
  private Node[] nodes;
  private Edge[] edges;
  private int width;
  private int height;
  
  public LayoutCalculator(
    Node[] nodes,
    Edge[] edges,
    int width,
    int height
  ){
    this.nodes = nodes;
    this.edges = edges;
    this.width = width;
    this.height = height;
    nnodes = nodes.length;
    nedges = edges.length;
  }

  public void layout(
    int iterates, 
    int movementLimit, 
    int gravity,
    int repulsionMultiplier
  ){
    //System.out.println("---------------");
    //System.out.println("--ATTRACTION--");
    for(int count=0; count<iterates; count++){
      //System.out.println("Count: "+count);
      
      //
      //Attraction
      for (int i = 0 ; i < nedges ; i++) {
        Edge e = edges[i];
        double vx = nodes[e.to].x - nodes[e.from].x;
        double vy = nodes[e.to].y - nodes[e.from].y;
        double len = Math.sqrt(vx * vx + vy * vy);
        double f = (edges[i].len - len) / (len * 3) ;
        double dx = f * vx;
        double dy = f * vy;
        
        //System.out.println("f: "+f+"\t"+dx+"\t"+dy);

        nodes[e.to].dx += dx;
        nodes[e.to].dy += dy;
        nodes[e.from].dx += -dx;
        nodes[e.from].dy += -dy;
        
      }

      //System.out.println("---------------");
      //System.out.println("---REPULSION--");
      
      //
      //Gravity
      if(gravity != 0){
        for (int i = 0 ; i < nnodes ; i++) {
          nodes[i].dy = nodes[i].dy - gravity;
        }
      }
      
      //
      //Repulsion
      for (int i = 0 ; i < nnodes ; i++) {
        Node n1 = nodes[i];
        double dx = 0;
        double dy = 0;

        for (int j = 0 ; j < nnodes ; j++) {
          if (i == j) {
            continue;
          }

          Node n2 = nodes[j];
          double vx = n1.x - n2.x;
          double vy = n1.y - n2.y;
          double len = vx * vx + vy * vy;
          if (len == 0) {
              dx += Math.random();
              dy += Math.random();
          } else if (len < width*height) {
              dx += vx / len;
              dy += vy / len;
          }
        }
        
        double dlen = dx * dx + dy * dy;
        if (dlen > 0) {
          dlen = Math.sqrt(dlen) / 2;
          n1.dx += repulsionMultiplier * dx / dlen;
          n1.dy += repulsionMultiplier * dy / dlen;
        }
        
        //System.out.println(n1.lbl+"\t"+n1.dx+"\t"+n1.dy+"\tlen:"+dlen);
        
      }

      //System.out.println("---------------");
      //System.out.println("--LIMITING----");
      //System.out.println("movement limit:"+movementLimit);

      Dimension d = new Dimension(width, height);
      for (int i = 0 ; i < nnodes ; i++) {
        Node n = nodes[i];
        //System.out.println("Before limit: "+n.lbl+"\t"+n.x+",\t "+n.y+" --\t"+n.dx+", \t"+n.dy);
        if (!n.fixed) {
          //System.out.println("x-Min of limit: "+Math.min(movementLimit, n.dx)+" Max: "+Math.max(-movementLimit, Math.min(movementLimit, n.dx)));
          n.x += Math.max(-movementLimit, Math.min(movementLimit, n.dx));
          n.y += Math.max(-movementLimit, Math.min(movementLimit, n.dy));
          if (n.x < 0) {
            if(n.dx < 0){
              n.x = -1*n.dx;
            }else{
              n.x = n.dx;
            }
          } else if (n.x > d.width) {
              n.x = d.width;
          }

          if (n.y < 0) {
            //System.out.println("neg ny");
            if(n.dy < 0){
              n.y = -1*n.dy;
            }else{
              n.y = n.dy;
            }
          } else if (n.y > d.height) {
              n.y = d.height;
          }
        }
        
        //System.out.println("\t"+n.lbl+":\t"+n.x+",\t"+n.y+" :: \t"+n.dx+",\t"+n.dy);
        n.dx /= 2;
        n.dy /= 2;
        //System.out.println("--------- "+n.lbl+" "+n.x+" , "+n.y+" -------- "+n.dx+" , "+n.dy);
      }
    }
  }
  
  public Node[] getNodes(){
    return nodes;
  }
}
