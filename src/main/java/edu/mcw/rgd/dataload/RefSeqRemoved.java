package edu.mcw.rgd.dataload;

import edu.mcw.rgd.process.FileDownloader;

/**
 * Created by IntelliJ IDEA.
 * User: mtutaj
 * Date: 10/3/14
 * Time: 3:57 PM
 * <p>
 * compare list of removed RefSeq acc ids against RGD and generate a report
 * </p>
 */
public class RefSeqRemoved {
    private String version;
    private String dataDir;
    private String dataDirAr;
    private String fileMask;
    private String reportFile;

    public void run() throws Exception {

        System.out.println(getVersion());

        FileDownloader downloader = new FileDownloader();
        downloader.setExternalFile("ftp://ftp.ncbi.nih.gov/refseq/removed");
        String[] files = downloader.listFiles();
        for( String file: files ) {
            System.out.println(file);
        }
    }

    public void setVersion(String version) {
        this.version = version;
    }

    public String getVersion() {
        return version;
    }

    public void setDataDir(String dataDir) {
        this.dataDir = dataDir;
    }

    public String getDataDir() {
        return dataDir;
    }

    public void setDataDirAr(String dataDirAr) {
        this.dataDirAr = dataDirAr;
    }

    public String getDataDirAr() {
        return dataDirAr;
    }

    public void setFileMask(String fileMask) {
        this.fileMask = fileMask;
    }

    public String getFileMask() {
        return fileMask;
    }

    public void setReportFile(String reportFile) {
        this.reportFile = reportFile;
    }

    public String getReportFile() {
        return reportFile;
    }
}
