import pandas as pd

def save_excel(df: pd.DataFrame, outname: str = "final.table", skip_adjust: list = []) -> None:
    outname: str = f"{outname}.xlsx" if "xlsx" not in outname else f"{outname}"
    sheetname: str = "Sheet1"

    writer = pd.ExcelWriter(outname, engine='xlsxwriter')
    df.to_excel(writer, index=False, sheet_name=sheetname, na_rep='NA')
    worksheet = writer.sheets[sheetname]  # pull worksheet object

    # Adjust columnss
    for idx, col in enumerate(df):  # loop through all columns
        series = df[col]
        max_len = max((
            series.astype(str).map(len).max(),  # len of largest item
            len(str(series.name))  # len of column name/header
            )) + 4  # adding a little extra space

        if col not in skip_adjust:
            worksheet.set_column(idx, idx, max_len)  # set column width
    
    writer.close()