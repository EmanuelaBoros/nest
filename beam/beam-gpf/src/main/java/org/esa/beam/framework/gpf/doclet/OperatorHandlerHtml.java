package org.esa.beam.framework.gpf.doclet;

import com.sun.javadoc.RootDoc;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

// todo - use template engine, e.g. apache velocity (nf)
public class OperatorHandlerHtml implements OperatorHandler {

    ArrayList<OperatorDesc> operatorDescs;
    private File baseDir;

    public OperatorHandlerHtml() {
        operatorDescs = new ArrayList<OperatorDesc>();
    }

    @Override
    public void start(RootDoc root) throws IOException {
        baseDir = new File("help");   // todo - use main arguments (nf)
        if (!baseDir.isDirectory()) {
            if (!baseDir.mkdirs()) {
                throw new IOException("Failed to create base directory " + baseDir);
            }
        }
    }

    @Override
    public void stop(RootDoc root) throws IOException {
        File indexFile = new File(baseDir, "index.html");
        PrintWriter writer = new PrintWriter(new FileWriter(indexFile));
        try {
            writeIndex(writer);
        } finally {
            writer.close();
        }
    }

    @Override
    public void processOperator(OperatorDesc operatorDesc) throws IOException {
        File file = getOperatorPageFile(operatorDesc);
        PrintWriter writer = new PrintWriter(new FileWriter(file));
        try {
            writeOperatorPage(operatorDesc, writer);
            operatorDescs.add(operatorDesc);
        } finally {
            writer.close();
        }
    }

    private void writeIndex(PrintWriter writer) throws IOException {
        writer.println("<!DOCTYPE HTML PUBLIC " +
                "\"-//W3C//DTD HTML 4.01//EN\" " +
                "\"http://www.w3.org/TR/html4/strict.dtd\">");
        writer.println("<html>");
        writer.println("<head><title>Operator Index</title></head>");
        writer.println("<body>");

        writer.println("<h1>Operator Index</h1>");

        writer.println("<table>");
        for (OperatorDesc operatorDesc : operatorDescs) {
            writer.println("  <tr>");
            writer.println("    <td><b><code><a href=\"" + getOperatorPageName(operatorDesc) + "\">" + operatorDesc.getName() + "</a></code></b></td>");
            writer.println("    <td>" + operatorDesc.getShortDescription() + "</td>");
            writer.println("  </tr>");
        }
        writer.println("</table>");
        writer.println("</body>");
        writer.println("</html>");
    }

    private void writeOperatorPage(OperatorDesc operatorDesc, PrintWriter writer) {
        writer.println("<!DOCTYPE HTML PUBLIC " +
                "\"-//W3C//DTD HTML 4.01//EN\" " +
                "\"http://www.w3.org/TR/html4/strict.dtd\">");
        writer.println("<html>");
        writer.println("<head><title>" + operatorDesc.getName() + "</title></head>");
        writer.println("<body>");

        writer.println("<h1>" + operatorDesc.getName() + " Operator</h1>");

        writer.println("<h2>Overview</h2>");
        writer.println("<table>");
        writer.println("  <tr><td><b>Name:</b></td><td><code>" + operatorDesc.getName() + "</code></td></tr>");
        writer.println("  <tr><td><b>Full name:</b></td><td><code>" + operatorDesc.getType().getName() + "</code></td></tr>");
        writer.println("  <tr><td><b>Purpose:</b></td><td>" + operatorDesc.getShortDescription() + "</td></tr>");
        writer.println("  <tr><td><b>Version:</b></td><td>" + operatorDesc.getVersion() + "</td></tr>");
        writer.println("</table>");

        writer.println("<h2>Description</h2>");
        String description = operatorDesc.getLongDescription();
        if (!description.isEmpty()) {
            writer.println(description);
        } else {
            writer.println("<i>No description available.</i>");
        }

        writer.println("<h2>Sources</h2>");
        SourceProductDesc[] sourceProducts = operatorDesc.getSourceProducts();
        if (sourceProducts.length > 0) {
            writer.println("<table>");
            writer.println("<tr>");
            writer.println("  <th>Name</th>");
            writer.println("  <th>Description</th>");
            writer.println("</tr>");
            for (SourceProductDesc sourceProduct : sourceProducts) {
                writer.println("<tr>");
                writer.println("  <td><code>" + sourceProduct.getName() + "</code></td>");
                writer.println("  <td>" + getFullDescription(sourceProduct) + "</td>");
                writer.println("</tr>");
            }
            writer.println("</table>");
        } else {
            writer.println("<p><i>This operator does not have any sources.</i></p>");
        }

        writer.println("<h2>Parameters</h2>");
        ParameterDesc[] parameterDescs = operatorDesc.getParameters();
        if (parameterDescs.length > 0) {
            writer.println("<table>");
            writer.println("<tr>");
            writer.println("  <th>Name</th>");
            writer.println("  <th>Type</th>");
            writer.println("  <th>Description</th>");
            writer.println("</tr>");
            for (ParameterDesc parameterDesc : parameterDescs) {
                writer.println("<tr>");
                writer.println("  <td><code>" + parameterDesc.getName() + "</code></td>");
                writer.println("  <td><code>" + parameterDesc.getField().getType().getSimpleName() + "</code></td>");
                writer.println("  <td>" + getFullDescription(parameterDesc) + "</td>");
                writer.println("</tr>");
            }
            writer.println("</table>");
        } else {
            writer.println("<p><i>This operator does not have any parameters.</i></p>");
        }

        writer.println("<h2>Usage</h2>");
        writer.println("<p><i>TODO</i></p>");

        writer.println("</body>");
        writer.println("</html>");
    }

    private String getFullDescription(ElementDesc elementDesc) {
        String shortDescription = elementDesc.getShortDescription();
        String longDescription = elementDesc.getLongDescription();
        if (shortDescription.isEmpty()) {
            return longDescription;
        }
        if (longDescription.isEmpty()) {
            return shortDescription;
        }
        return shortDescription + "<br/>" + longDescription;
    }

    private File getOperatorPageFile(OperatorDesc operatorDesc) {
        return new File(baseDir, getOperatorPageName(operatorDesc));
    }

    private String getOperatorPageName(OperatorDesc operatorDesc) {
        return operatorDesc.getType().getName().replace('.', '_') + ".html";
    }

}