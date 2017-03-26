VERSION 1.0 CLASS
BEGIN
  MultiUse = -1  'True
END
Attribute VB_Name = "ThisWorkbook"
Attribute VB_GlobalNameSpace = False
Attribute VB_Creatable = False
Attribute VB_PredeclaredId = True
Attribute VB_Exposed = True
' Set string comparisons to be case insensitive
Option Compare Text

Sub Main()
    ' Turn off screen updating to speed performance
    Application.ScreenUpdating = False

    ' Set up handle for macro file
    Dim macroFile As Workbook
    Set macroFile = ThisWorkbook
    Dim metadataSheet As Worksheet
    Set metadataSheet = macroFile.Sheets(1)
    
    ' Specify input files/folders
    Dim inputFiles As Variant
    Set inputFiles = PromptForFiles(True)
    
    ' Create new output file
    Dim outputFile As Workbook
    Set outputFile = CreateNewWorkbook()
    
    ' Iterate each input file
    For Each inputFileName In inputFiles
    
        ' Open the input file
        Dim inputFile As Workbook
        Set inputFile = Workbooks.Open(inputFileName)
        
        ' Iterate each label pairing
        With metadataSheet
            numPairs = .Cells(1, 1).End(xlDown).Row
            For pairIdx = 2 To numPairs Step 1
                
                ' Copy relevant data to output file
                Dim label1 As String, label2 As String, metric As String
                label1 = .Cells(pairIdx, 1)
                label2 = .Cells(pairIdx, 2)
                metric = .Cells(pairIdx, 3)
                CopyData inputFile, outputFile, label1, label2, metric
            Next pairIdx
        End With
        
        ' Close the input file
        Application.DisplayAlerts = False
        inputFile.Close
        Application.DisplayAlerts = True
    Next inputFileName
    
    ' Delete Sheet1
    Application.DisplayAlerts = False
    outputFile.Sheets("Sheet1").Delete
    Application.DisplayAlerts = True
    
    ' Add final formulas
    With metadataSheet
        numPairs = .Cells(1, 1).End(xlDown).Row
        For pairIdx = 2 To numPairs Step 1
            label1 = .Cells(pairIdx, 1)
            label2 = .Cells(pairIdx, 2)
            metric = .Cells(pairIdx, 3)
            sheetLabel1 = Split(label1, ";")(0)
            sheetLabel2 = Split(label2, ";")(0)
            sheetName = sheetLabel1 + " vs " + sheetLabel2
            Dim sheet As Worksheet
            Set sheet = outputFile.Sheets(sheetName)
            PopulateAnalysisFormulas sheet, label1, label2, metric
        Next pairIdx
    End With
    
    Application.ScreenUpdating = True
    MsgBox "Success"
End Sub

Function PromptForFiles(isMultiple As Boolean) As Variant

    ' Init
    Dim dialog As FileDialog
    Set dialog = Application.FileDialog(msoFileDialogFilePicker)
            
    With dialog
        ' Allow multiple files to be selected, if desired
        .AllowMultiSelect = isMultiple
    
        ' Only show excel files
        .Filters.Add "Excel", "*.xls; *.xlsx", 1
    
        ' Show the file picker dialog
        If .Show = -1 Then
            ' Return the files selected
            Set PromptForFiles = .SelectedItems
        End If
    End With
    
    ' Cleanup
    Set dialog = Nothing
End Function

Function CreateNewWorkbook() As Variant

    ' Setup the new workbook
    Dim newBook As Workbook
    Set newBook = Workbooks.Add
    With newBook
    
        ' Delete all but one worksheet
        Application.DisplayAlerts = False
        For i = .Sheets.Count To 2 Step -1
            .Sheets(i).Delete
        Next i
        Application.DisplayAlerts = True
        
        ' Activate the workbook
        .Activate
    End With
    
    ' Return the new workbook
    Set CreateNewWorkbook = newBook
End Function

Function ReadMacroMetadata() As String

    ' Select the macro workbook
    Dim macroWorkbook As Workbook
    Set macroWorkbook = Workbooks("Macro.xlsm")
    With macroWorkbook
        
        ' Select the worksheet
        With macroWorkbook.Sheets(1)
        
            ' Grab the whole box of data, minus headers
            Dim myRange As Variant
            ' myRange = .Range(.Cells(2, 1), .Cells(.Cells(2, 1).End(xlDown).Row, 3)).Address
            topLeft = .Cells(2, 1).Address
            bottomRight = .Cells(.Range(topLeft).End(xlDown).Row, 3).Address
            ReadMacroMetadata = .Range(topLeft, bottomRight).Address
        End With
    End With
    
    ' Cleanup
    Set macroWorkbook = Nothing
End Function

Sub CopyData(inputBook As Workbook, outputBook As Workbook, label1 As String, label2 As String, metric As String)

    ' Set up the output worksheet
    Dim outputSheet As Worksheet
    Set outputSheet = SetUpOutputSheet(inputBook, outputBook, label1, label2)
    
    ' Set up row pointers
    label1InputRow = 0
    label2InputRow = 0
    lastRow = outputSheet.Cells(outputSheet.Rows.Count, 2).End(xlUp).Row
    label1OutputRow = lastRow + 1
    label2OutputRow = lastRow + 2
    
    ' Get rows to copy for label1
    With inputBook.Sheets(1)
        maxRow = .Cells(2, 1).End(xlDown).Row
        
        ' Identify the row in input that matches label1
        For currentRow = 2 To maxRow Step 1
            If MatchLabels(.Cells(currentRow, 2).Value, label1) Then
            
                ' Match found - grab the row number and break
                label1InputRow = currentRow
                Exit For
            End If
        Next currentRow
        
        ' Get the label2 row as well
        For currentRow = 2 To maxRow Step 1
            If MatchLabels(.Cells(currentRow, 2).Value, label2) Then
            
                label2InputRow = currentRow
                Exit For
            End If
        Next currentRow
    End With
    
    ' Make sure both labels exist in input file before copying
    If label1InputRow <> 0 And label2InputRow <> 0 Then
        
        ' Copy the rows from input to output
        inputBook.Sheets(1).Rows(label1InputRow).Copy
        outputSheet.Cells(label1OutputRow, 1).PasteSpecial
        inputBook.Sheets(1).Rows(label2InputRow).Copy
        outputSheet.Cells(label2OutputRow, 1).PasteSpecial
        
        ' Add a whitespace character in the next row to block it off
        outputSheet.Cells(label2OutputRow + 1, 2).Value = " "
        
        ' Add Ratio of MGI
        metricColumn = GetMetricColumn(outputSheet, metric)
        ratioColumn = outputSheet.Cells(1, 1).End(xlToRight).column + 2
        val1 = outputSheet.Cells(label1OutputRow, metricColumn).Address
        val2 = outputSheet.Cells(label2OutputRow, metricColumn).Address
        outputSheet.Cells(label1OutputRow, ratioColumn).Formula = "=" + val1 + "/" + val2
        
        ' Add header for ratio column, if it doesn't exist
        If label1OutputRow = 2 Then
            outputSheet.Cells(1, ratioColumn).Value = "Ratio of MGI"
        End If
    End If
End Sub

Function MatchLabels(inputLabel As String, desiredLabels As String) As Boolean

    ' Parse the input label
    parsedInputLabel = Split(inputLabel, ":")(1)
    
    ' Iterate each label option in desired labels
    For Each labelOption In Split(desiredLabels, ";")
    
        ' Compare the labels and break early if there is a match
        If parsedInputLabel = labelOption Then
            MatchLabels = True
            Exit Function
        End If
    Next labelOption
    
    MatchLabels = False
End Function

Function SetUpOutputSheet(inputBook As Workbook, outputBook As Workbook, label1, label2)

    ' Generate worksheet name from labels
    sheetLabel1 = Split(label1, ";")(0)
    sheetLabel2 = Split(label2, ";")(0)
    sheetName = sheetLabel1 + " vs " + sheetLabel2
    
    ' Check if worksheet already exists
    outputBook.Activate
    If Not Evaluate("ISREF('" & sheetName & "'!A1)") Then
        
        ' Worksheet does not exist - create it
        Dim newSheet As Worksheet
        With outputBook
            Set newSheet = .Sheets.Add(After:=.Sheets(.Sheets.Count))
            newSheet.Name = sheetName
        End With
        
        ' Copy the headers
        inputBook.Sheets(1).Rows(1).Copy
        newSheet.Range("A1").PasteSpecial
        newSheet.Rows(1).Font.Bold = True
        
        ' Return the newly created worksheet
        Set SetUpOutputSheet = newSheet
    Else
        
        ' Worksheet already exists - return that one
        Set SetUpOutputSheet = outputBook.Sheets(sheetName)
    End If
End Function

Function GetMetricColumn(sheet As Worksheet, metricName As String) As Integer

    lastColumn = sheet.Cells(1, 1).End(xlToRight).column
    For col = 2 To lastColumn Step 1
        If sheet.Cells(1, col).Value = metricName Then
            GetMetricColumn = col
            Exit Function
        End If
    Next col
    
    GetMetricColumn = 0
End Function

Sub PopulateAnalysisFormulas(sheet As Worksheet, label1 As String, label2 As String, metricName As String)

    ' Average MGI for label1
    label1Column = sheet.Cells(1, 2).End(xlToRight).column + 4
    label1Cells = GetListOfCells(sheet, label1, GetMetricColumn(sheet, metricName))
    label1AvgHeader = "Average " + Split(label1, ";")(0) + " MGI"
    label1AvgFunction = "=AVERAGE(" + label1Cells + ")"
    sheet.Cells(1, label1Column).Value = label1AvgHeader
    sheet.Cells(1, label1Column).Font.Bold = True
    sheet.Cells(2, label1Column).Value = label1AvgFunction
    
    ' Std Dev MGI for label1
    label1StdHeader = "Std " + Split(label1, ";")(0) + " MGI"
    label1StdFunction = "=STDEV(" + label1Cells + ")"
    sheet.Cells(3, label1Column).Value = label1StdHeader
    sheet.Cells(3, label1Column).Font.Bold = True
    sheet.Cells(4, label1Column).Value = label1StdFunction
    
    ' Average MGI for label2
    label2Column = label1Column + 2
    label2Cells = GetListOfCells(sheet, label2, GetMetricColumn(sheet, metricName))
    label2AvgHeader = "Average " + Split(label2, ";")(0) + "MGI"
    label2AvgFunction = "=AVERAGE(" + label2Cells + ")"
    sheet.Cells(1, label2Column).Value = label2AvgHeader
    sheet.Cells(1, label2Column).Font.Bold = True
    sheet.Cells(2, label2Column).Value = label2AvgFunction
    
    ' Std Dev MGI for label2
    label2StdHeader = "Std " + Split(label2, ";")(0) + " MGI"
    label2StdFunction = "=STDEV(" + label2Cells + ")"
    sheet.Cells(3, label2Column).Value = label2StdHeader
    sheet.Cells(3, label2Column).Font.Bold = True
    sheet.Cells(4, label2Column).Value = label2StdFunction
    
    ' Ratio of Average MGI
    ratioColumn = label2Column + 2
    ratioHeader = "Ratio of Average MGI"
    ratioFunction = "=" & sheet.Cells(2, label1Column).Address & "/" & sheet.Cells(2, label2Column).Address
    sheet.Cells(1, ratioColumn).Value = ratioHeader
    sheet.Cells(1, ratioColumn).Font.Bold = True
    sheet.Cells(2, ratioColumn).Value = ratioFunction
    
    ' Average of the Ratios
    ratioAvgHeader = "Average of Individual Ratios"
    ratioAvgCells = GetListOfCells(sheet, label1, sheet.Cells(1, 2).End(xlToRight).column + 2)
    ratioAvgFunction = "=AVERAGE(" & ratioAvgCells & ")"
    sheet.Cells(3, ratioColumn).Value = ratioAvgHeader
    sheet.Cells(3, ratioColumn).Font.Bold = True
    sheet.Cells(4, ratioColumn).Value = ratioAvgFunction
End Sub

Function GetListOfCells(sheet As Worksheet, label As String, column As Integer) As String

    ' Build a comma-delimited list of cells containing label
    Dim runningList As String
    
    ' Iterate each row (again, ugh)
    lastRow = sheet.Cells(sheet.Rows.Count, 2).End(xlUp).Row
    For currentRow = 2 To lastRow Step 1
        If Trim(sheet.Cells(currentRow, 2).Value + vbNullString) <> vbNullString Then
            If MatchLabels(sheet.Cells(currentRow, 2).Value, label) Then
                ' Found a match - add it to the running list
                If runningList = "" Then
                    runningList = runningList + sheet.Cells(currentRow, column).Address
                Else
                    runningList = runningList + "," + sheet.Cells(currentRow, column).Address
                End If
            End If
        End If
    Next currentRow
    
    GetListOfCells = runningList
End Function

Function GetMetadata(sheet As Worksheet, label1 As String, label2 As String) As String

    lastRow = sheet.Cells(1, 1).End(xlDown).Row
    For currentRow = 2 To lastRow Step 1
        If MatchLabels("x:" + label1, sheet.Cells(currentRow, 1)) And MatchLabels("x:" + label2, sheet.Cells(currentRow, 2)) Then
            GetMetadata = sheet.Cells(currentRow, 3)
            Exit Function
        End If
    Next currentRow
End Function
